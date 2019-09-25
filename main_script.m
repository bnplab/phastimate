%% TODO

% [X] true amplitude, trial by trial data (adapt phastimate_truephase to also output amplitude)
% [X] visualize equivalent filters, optimize filter object selection
% [X] estimate phase error with Zrenner2018 parameters
% [ ] automated artifact rejection, compare results for data with rejected artifacts
% [ ] include graphical output option in phastimate (for debugging and explanation, including as figure)
% [ ] move filter object creation to separate phastimate_create_filter_objects.m
% [ ] refactor Step 4
% [ ] include more subjects (?)

%%

% Known issues/Developtment TODO:
% - peak frequency has a resolution of 0.5Hz, could be improved
% - change create_epochs to create predetermined number of epochs, currently NUMEPOCHS has to match this exactly
% - when creating subplot figures, subplot should refer to current figure to avoid change of current figure by user moving plot targets
% - there is no protection against designing filters that require longer windows than the available epoch length (filter order depends on peak frequency, range is set in PEAK_FREQUENCY_INTERVAL constant)

%% This is the script that generates the data for the figures and demonstrates the code

% Note:
%  - the script uses estimate_SNR.m to fit 1/f noise and determine SNR
%  - circular variance is used as the measure of dispersion (note that circular variance ranges between 0 and 1 and is equal to 1-PLV)
%  - all circular data is in radians
%  - peak frequency is used to individualize filters in phastimate (0.5 Hz resolution)
%  - epoch creation happens when needed and needs to be identical in different parts of the script
%  - Signal Processing and Global Optimization Toolboxes are required
%  - This script has been tested with Matlab 2017b

%% preliminaries

% check for toolboxes
assert(~isempty(which('designfilt')), 'filter design function designfilt.m not found, is the Signal Processing Toolbox installed?')
assert(~isempty(which('range')), 'statistical function range.m not found, is the Statistics and Machine Learning Toolbox installed?')
assert(~isempty(which('ga')), 'genetic algorithm function ga.m not found, is the Global Optimization Toolbox installed?')

% clear variables, close windows, reset paths
clear all; close all; path(pathdef); clc

% switch to current directory and add relative path to phastimate toolbox
cd(fileparts(getfield(matlab.desktop.editor.getActive, 'Filename')))
addpath('../phastimate_code/')

% circular statistics functions (simplified from circstat toolbox)
ang_mean = @(x) angle(mean(exp(1i*x)));
ang_diff = @(x, y) angle(exp(1i*x)./exp(1i*y));
ang_var = @(x) 1-abs(mean(exp(1i*x)));
%ang_var2dev = @(v) sqrt(2*v); % circstat preferred formula uses angular deviation (bounded from 0 to sqrt(2)) which is sqrt(2*(1-r))
ang_var2dev = @(v) sqrt(-2*log(1-v)); % formula for circular standard deviation is sqrt(-2*ln(r))

%% set constants

% filter design method for phastimate (order and peak frequency is variable)
design_phastimate_filter = @(ord, freq, fs) designfilt('bandpassfir', 'FilterOrder', ord, 'CutoffFrequency1', freq-1, 'CutoffFrequency2', freq+1, 'SampleRate', fs, 'DesignMethod', 'window');

NUM_EPOCHS = 497;
HILBERTWIN = 128; % this is an appropriate window for alpha at 1000 Hz
PEAK_FREQUENCY_INTERVAL = [8 14];

%% load resting sate data into a master data table 'T'

data = load(fullfile('..', 'data', 'murhythmdataset.mat')); 

T = table('RowNames', data.subject_ids);
T.data = data.data;
T.fs = data.fs * ones(height(T),1);

clear('data');


%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine alpha peak frequency and signal to noise ratio
% Note:
% - data is not cleaned before this step
% - figure is created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T.peak_frequency = nan(height(T),1);
T.peak_SNR = nan(height(T),1);

figures_snr = {};
figures_snr{1} = figure;

subplot_index = 1;

for row_index = 1:height(T)
    subject = T(row_index,:).Row{1};
    
    epochs = create_epochs_overlapping(T(row_index,:).data, T(row_index,:).fs); % from continuous data

    % estimate SNR and plot
    ax = subplot(10,5,subplot_index);
    [peak_frequency, peak_SNR] = estimate_SNR(epochs, T(row_index,:).fs, PEAK_FREQUENCY_INTERVAL, ax);
     
    % save data
    if ~isempty(peak_frequency)
        T(row_index,:).peak_frequency = peak_frequency;
        T(row_index,:).peak_SNR = peak_SNR;
    end
    
    % beautify axes
    legend off
    title(ax, {subject; ax.Title.String}, 'Interpreter', 'none')
    set(ax, 'XTick', [3, 5, 8, 14, 30])
    if (subplot_index < 46 && row_index <= height(T)-5), xlabel(ax,''); end % unless bottom row
    if (mod(subplot_index,5) ~= 1), ylabel(ax,''); end % cunless first column
    
    subplot_index = subplot_index + 1;
    drawnow
    
    if subplot_index > 50 % new page?
        figures_snr{numel(figures_snr)+1} = figure;
        subplot_index = 1;
    end
end

% save figures as PDF
for i = 1:numel(figures_snr)
    set(figures_snr{i},'Renderer','Painters') %export vectorized
    set(figures_snr{i},'PaperPositionMode','manual','PaperType','a2','PaperOrientation','portrait','PaperUnits','normalized','PaperPosition',[0 0 1 1])
    print(figures_snr{i},['figure_snr (page ' num2str(i) ')'], '-dpdf', '-r0')
end

% exlude subjects without a peak or with negative peak SNR

fprintf('\nRemoving %i/%i entries where no peak could be found [ %s]... ', sum(isnan(T.peak_frequency)), height(T), sprintf('%s ', string(T.Row(isnan(T.peak_frequency)))))
T(isnan(T.peak_frequency), :) = [];
fprintf('done')
fprintf('\nRemoving %i/%i entries with a negative SNR [ %s]... ', sum(T.peak_SNR < 0), height(T), sprintf('%s ', string(T.Row(T.peak_SNR < 0))))
T(T.peak_SNR < 0, :) = [];
fprintf('done')
fprintf('\nThe number of remaining datasets is %i\n', height(T))

clear('subplot_index', 'row_index', 'subject', 'ax', 'epochs', 'peak_frequency', 'peak_SNR', 'i')


%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine true phase and amplitude, as well as variance of "true" phase
% - epochs are recreated from the data when needed to save memory space
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add columns for epoch-by-epoch data
T.epochs_truephase_mean = nan(height(T), NUM_EPOCHS);
T.epochs_truephase_ang_var = nan(height(T), NUM_EPOCHS);
T.epochs_trueamp_mean = nan(height(T), NUM_EPOCHS);
T.epochs_trueamp_var = nan(height(T), NUM_EPOCHS);

figures_truephase = {};
figures_truephase{1} = figure;

subplot_index = 1;

for row_index = 1:height(T)
    subject = T(row_index,:).Row{1};
    
    epochs = create_epochs_overlapping(T(row_index,:).data, T(row_index,:).fs); % from continuous data

    peak_frequency = T(row_index,:).peak_frequency;
    
    % set-up equivalent filter objects for given peak frequency
    filter_objects = {};
    fs = T(row_index,:).fs;
    for ord = [2 3 4 5] % FIR - windowed sinc
        filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'CutoffFrequency1', peak_frequency-1, 'CutoffFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'window')};
    end
    for ord = [3 4 5] % FIR - least squares (equiripple is similar)
        filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'StopbandFrequency1', peak_frequency-4, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+4, 'SampleRate', fs, 'DesignMethod', 'ls')};
    end
    for ord = [4 8 12] % IIR - butterworth
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'HalfPowerFrequency1', peak_frequency-1, 'HalfPowerFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'butter')};
    end
    for ord = [4 6 8] % IIR - chebychev I
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'PassbandFrequency1', peak_frequency-1.5, 'PassbandFrequency2', peak_frequency+1.5, 'PassbandRipple', 0.5, 'SampleRate', fs, 'DesignMethod', 'cheby1')};
    end
    for attenuation = [10 20] % IIR - elliptic
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'StopbandFrequency1', peak_frequency-2, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+2, 'StopbandAttenuation1', attenuation, 'PassbandRipple', 0.5, 'StopbandAttenuation2', attenuation, 'SampleRate', fs, 'DesignMethod', 'ellip', 'MatchExactly', 'passband')};
    end    
    
    [truephase_mean, truephase_variance, trueamp_mean, trueamp_variance] = phastimate_truephase(epochs, filter_objects);
    
    T(row_index,:).epochs_truephase_mean = truephase_mean;
    T(row_index,:).epochs_truephase_ang_var = truephase_variance;
    
    T(row_index,:).epochs_trueamp_mean = trueamp_mean;
    T(row_index,:).epochs_trueamp_var = trueamp_variance;
    
    ax = subplot(10,5,subplot_index); hold on
    histogram(ax, rad2deg(ang_var2dev(T(row_index,:).epochs_truephase_ang_var)), 'BinWidth', 1, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
    plot(ax, [1 1] .* rad2deg(ang_var2dev(quantile(T(row_index,:).epochs_truephase_ang_var, 0.5))), [0 0.5], 'LineWidth', 2, 'Color', 'red')
    xlim(ax, [0 20]); ylim(ax, [0 1]);
    if (subplot_index > 45 || row_index > height(T)-5), xlabel(ax,'circular deviation (deg)'); end % bottom row
    if (mod(subplot_index,5) == 1), ylabel(ax,'cumulative probability'); end % first column
    title(subject, 'Interpreter', 'none')
    
    subplot_index = subplot_index + 1;
    drawnow
    
    if subplot_index > 50
        figures_truephase{numel(figures_truephase)+1} = figure;
        subplot_index = 1;
    end

end

% save figures as PDF
for i = 1:numel(figures_truephase);
    set(figures_truephase{i},'Renderer','Painters') %export vectorized
    set(figures_truephase{i},'PaperPositionMode','manual','PaperType','a2','PaperOrientation','portrait','PaperUnits','normalized','PaperPosition',[0 0 1 1])
    print(figures_truephase{i},['figure_truephase (page ' num2str(i) ')'], '-dpdf', '-r0')
end

clear('subplot_index', 'row_index', 'subject', 'epochs', 'peak_frequency', 'filter_objects', 'fs', 'ord', 'truephase_mean', 'truephase_variance', 'trueamp_mean', 'trueamp_variance', 'ax', 'h', 'i')


%% Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine optimized phastimate parameters and resulting estimate
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add relevant columns
T.optim_window = nan(height(T), 1);
T.optim_filter_ord = nan(height(T), 1);
T.optim_edge = nan(height(T), 1);
T.optim_ar_ord = nan(height(T), 1);
T.optim_fval = nan(height(T), 1);

fprintf('\nNow running genetic algorithm to find optimized phastimate parameters...')
for row_index = 1:height(T) % iterate through entries
    subject = T(row_index,:).Row{1};
    fprintf('\nProcessing %s ... ', subject);
    assert(T(row_index,:).fs == 1000, 'default bounds for genetic optimization algorithm are set for detecting 8-14 Hz alpha assuming a sample rate of 1kHz');

    epochs = create_epochs_overlapping(T(row_index,:).data, T(row_index,:).fs); % from continuous data

    peak_frequency = T(row_index,:).peak_frequency;

    filter_order_range = 100:250;

    filter_objects_by_order = {}; %the index has to correspond to the order for the genetic algorithm
    for ord = filter_order_range
        filter_objects_by_order{ord} = design_phastimate_filter(ord, peak_frequency, T(row_index,:).fs);
    end
    
    bounds_filter_order = [filter_order_range(1) filter_order_range(end)];
    bounds_window = [400 750];
    bounds_edge = [30 120];
    bounds_ar_order = [5 60];

    % the includemask allows optimizing for a subset of epochs
    % it makes sense to exclude epochs that would also be excluded by the
    % real-time system, e.g. if artifacts are detected so as to not optimize
    % for noisy epochs that wouldn't result in a stimulus anyway
    
    % subselect according to truephase variance
    %includemask = T(row_index,:).epochs_truephase_angdev <= quantile(T(row_index,:).epochs_truephase_angdev, 0.5);

    % subselect according to true amplitude
    includemask = T(row_index,:).epochs_trueamp_mean >= quantile(T(row_index,:).epochs_trueamp_mean, 0.5);
    
    [optimal_parameters, ga_output] = phastimate_optimize(epochs(1:ceil(end/2),includemask), T(row_index,:).epochs_truephase_mean(includemask), filter_objects_by_order, bounds_filter_order, bounds_window, bounds_edge, bounds_ar_order, HILBERTWIN);

    % rerun phastimate with the optimized settings to confirm result
    D = design_phastimate_filter(optimal_parameters.filter_order, peak_frequency, T(row_index,:).fs);
    [estphase, estamp] = phastimate(epochs(((-optimal_parameters.window_length+1):0)+ceil(end/2),:), D, optimal_parameters.edge, optimal_parameters.ar_order, 128);
    
    % sanity check if the angular deviation matches the result of the optimization
    phases_error = ang_diff(T(row_index,:).epochs_truephase_mean, estphase);
    assert(abs(optimal_parameters.fval - ang_var(phases_error(includemask))) < 0.01, 'could not replicate result of optimization, were the same filters used?')
    
    %TODO: save output of ga, including number of generations etc.

    T(row_index,:).optim_window = optimal_parameters.window_length;
    T(row_index,:).optim_filter_ord = optimal_parameters.filter_order;
    T(row_index,:).optim_edge = optimal_parameters.edge;
    T(row_index,:).optim_ar_ord = optimal_parameters.ar_order;
    T(row_index,:).optim_fval = optimal_parameters.fval;

end
fprintf('\nDone.\n')

clear('row_index', 'subject', 'epochs', 'peak_frequency', 'filter_order_range', 'ord', 'bounds_ar_order', 'filter_objects_by_order', 'bounds_edge', 'bounds_filter_order', 'bounds_window', 'D', 'includemask', 'estamp', 'estphase', 'optimal_parameters', 'ga_output', 'phases_error')

%% save data

fprintf('\nSaving the results table...')
save('results_table', 'T', '-v7.3')
fprintf('\nDone.\n')

% NOTE: It's possible to load the data and continue the script below
T(ismember(T.Row, {'SCREEN_001', 'SCREEN_017', 'SCREEN2_085'}),:) = [];
T.epochs_truephase_ang_var = (T.epochs_truephase_angdev.^2)/2;


%% Step 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate phase with various parametersets and determine epoch-by-epoch
% metrics
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine default parameters (median of optimized values)
default_parameters = [];
default_parameters.window_length = ceil(median(T.optim_window));
default_parameters.filter_order = ceil(median(T.optim_filter_ord));
default_parameters.edge = ceil(median(T.optim_edge));
default_parameters.ar_order = ceil(median(T.optim_ar_ord));

%add relevant columns
T.epochs_estphase_default = nan(height(T), NUM_EPOCHS);
T.epochs_estamp_default = nan(height(T), NUM_EPOCHS);
T.error_ang_var_default = nan(height(T), 1);
T.error_ang_mean_default = nan(height(T), 1);
T.epochs_param_range_default = nan(height(T), NUM_EPOCHS);
T.epochs_param_stdev_default = nan(height(T), NUM_EPOCHS);

T.epochs_estphase_default_individualpassband = nan(height(T), NUM_EPOCHS);
T.epochs_estamp_default_individualpassband = nan(height(T), NUM_EPOCHS);
T.error_ang_var_default_individualpassband = nan(height(T), 1);
T.error_ang_mean_default_individualpassband = nan(height(T), 1);
T.epochs_param_range_default_individualpassband = nan(height(T), NUM_EPOCHS);
T.epochs_param_stdev_default_individualpassband = nan(height(T), NUM_EPOCHS);

T.epochs_estphase_optim = nan(height(T), NUM_EPOCHS);
T.epochs_estamp_optim = nan(height(T), NUM_EPOCHS);
T.error_ang_var_optim = nan(height(T), 1);
T.error_ang_mean_optim = nan(height(T), 1);
T.epochs_param_range_optim = nan(height(T), NUM_EPOCHS);
T.epochs_param_stdev_optim = nan(height(T), NUM_EPOCHS);

T.error_ang_var_zrenner2018 = nan(height(T), 1);
T.error_ang_mean_zrenner2018 = nan(height(T), 1);
T.epochs_estphase_zrenner2018 = nan(height(T), NUM_EPOCHS);
T.epochs_estamp_zrenner2018 = nan(height(T), NUM_EPOCHS);
T.epochs_param_range_zrenner2018 = nan(height(T), NUM_EPOCHS);
T.epochs_param_stdev_zrenner2018 = nan(height(T), NUM_EPOCHS);

T.error_ang_var_zrenner2018_individualpassband = nan(height(T), 1);
T.error_ang_mean_zrenner2018_individualpassband = nan(height(T), 1);
T.epochs_estphase_zrenner2018_individualpassband = nan(height(T), NUM_EPOCHS);
T.epochs_estamp_zrenner2018_individualpassband = nan(height(T), NUM_EPOCHS);
T.epochs_param_range_zrenner2018_individualpassband = nan(height(T), NUM_EPOCHS);
T.epochs_param_stdev_zrenner2018_individualpassband = nan(height(T), NUM_EPOCHS);

fprintf('\nNow running phastimate with different sets of parameters...')
for row_index = 1:height(T) % iterate through entries
    subject = T(row_index,:).Row{1};
    fprintf('\nProcessing %s ... ', subject);

    epochs = create_epochs_overlapping(T(row_index,:).data, T(row_index,:).fs); % from continuous data
    peak_frequency = T(row_index,:).peak_frequency;
    
    fprintf('peak frequency %4.1f Hz ... ', peak_frequency);
    
    % default parameters phastimate, fixed filter passband 8..13Hz
    epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);
    
    D = designfilt('bandpassfir', 'FilterOrder', default_parameters.filter_order, 'CutoffFrequency1', 8, 'CutoffFrequency2', 13, 'SampleRate', T(row_index,:).fs, 'DesignMethod', 'window');
    [estphase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, HILBERTWIN);
    
    T(row_index,:).epochs_estphase_default = estphase;
    T(row_index,:).epochs_estamp_default = estamp;    
    
    phases_error = ang_diff(T(row_index,:).epochs_truephase_mean, estphase);
    T(row_index,:).error_ang_var_default = ang_var(phases_error);
    T(row_index,:).error_ang_mean_default = ang_mean(phases_error);
    
    T(row_index,:).epochs_param_range_default = range(epochs(epochwindowmask,:));
    T(row_index,:).epochs_param_stdev_default = std(epochs(epochwindowmask,:));
    
    % default parameters phastimate, with individual filter passband
    epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);
    
    D = design_phastimate_filter(default_parameters.filter_order, peak_frequency, T(row_index,:).fs);
    [estphase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, HILBERTWIN);
    
    T(row_index,:).epochs_estphase_default_individualpassband = estphase;
    T(row_index,:).epochs_estamp_default_individualpassband = estamp;    
    
    phases_error = ang_diff(T(row_index,:).epochs_truephase_mean, estphase);
    T(row_index,:).error_ang_var_default_individualpassband = ang_var(phases_error);
    T(row_index,:).error_ang_mean_default_individualpassband = ang_mean(phases_error);
    
    T(row_index,:).epochs_param_range_default_individualpassband = range(epochs(epochwindowmask,:));
    T(row_index,:).epochs_param_stdev_default_individualpassband = std(epochs(epochwindowmask,:));
    
    % individually optimized phastimate
    epochwindowmask = ((-T(row_index,:).optim_window+1):0)+ceil(size(epochs,1)/2);
    
    D = design_phastimate_filter(T(row_index,:).optim_filter_ord, peak_frequency, T(row_index,:).fs);
    [estphase, estamp] = phastimate(epochs(epochwindowmask,:), D, T(row_index,:).optim_edge, T(row_index,:).optim_ar_ord, HILBERTWIN);

    T(row_index,:).epochs_estphase_optim = estphase;
    T(row_index,:).epochs_estamp_optim = estamp;
    
    phases_error = ang_diff(T(row_index,:).epochs_truephase_mean, estphase);
    T(row_index,:).error_ang_var_optim = ang_var(phases_error);
    T(row_index,:).error_ang_mean_optim = ang_mean(phases_error);
    
    T(row_index,:).epochs_param_range_optim = range(epochs(epochwindowmask,:));
    T(row_index,:).epochs_param_stdev_optim = std(epochs(epochwindowmask,:));

    % phastimates with parameters used in Zrenner2018 paper
    epochwindowmask = ((-500+1):0)+ceil(size(epochs,1)/2); %500ms window

    assert(T(row_index,:).fs == 1000, 'phastimate with Zrenner2018 parameters requires sample rate of 1000 Hz')
    D = designfilt('bandpassfir', 'FilterOrder', 128, 'CutoffFrequency1', 8, 'CutoffFrequency2', 13, 'SampleRate', T(row_index,:).fs, 'DesignMethod', 'window');
    [estphase, estamp] = phastimate(epochs(epochwindowmask,:), D, 64, 30, HILBERTWIN);

    T(row_index,:).epochs_estphase_zrenner2018 = estphase;
    T(row_index,:).epochs_estamp_zrenner2018 = estamp;
    
    phases_error = ang_diff(T(row_index,:).epochs_truephase_mean, estphase);
    T(row_index,:).error_ang_var_zrenner2018 = ang_var(phases_error);
    T(row_index,:).error_ang_mean_zrenner2018 = ang_mean(phases_error);
    
    T(row_index,:).epochs_param_range_zrenner2018 = range(epochs(epochwindowmask,:));
    T(row_index,:).epochs_param_stdev_zrenner2018 = std(epochs(epochwindowmask,:));

    % phastimates with parameters used in Zrenner2018 paper, but with indiviudal peak frequency filter passband
    epochwindowmask = ((-500+1):0)+ceil(size(epochs,1)/2); %500ms window

    assert(T(row_index,:).fs == 1000, 'phastimate with Zrenner2018 parameters requires sample rate of 1000 Hz')
    D = design_phastimate_filter(128, peak_frequency, T(row_index,:).fs);
    [estphase, estamp] = phastimate(epochs(epochwindowmask,:), D, 64, 30, HILBERTWIN);

    T(row_index,:).epochs_estphase_zrenner2018_individualpassband = estphase;
    T(row_index,:).epochs_estamp_zrenner2018_individualpassband = estamp;
    
    phases_error = ang_diff(T(row_index,:).epochs_truephase_mean, estphase);
    T(row_index,:).error_ang_var_zrenner2018_individualpassband = ang_var(phases_error);
    T(row_index,:).error_ang_mean_zrenner2018_individualpassband = ang_mean(phases_error);
    
    T(row_index,:).epochs_param_range_zrenner2018_individualpassband = range(epochs(epochwindowmask,:));
    T(row_index,:).epochs_param_stdev_zrenner2018_individualpassband = std(epochs(epochwindowmask,:));

end
fprintf('\nDone.\n')

clear('row_index', 'subject', 'epochs', 'epochwindowmask', 'peak_frequency', 'D', 'estamp', 'estphase', 'phases_error')


%% Step 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization and Analysis
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Characteristics of signal spectra

figure

subplot(1,2,1)
histogram(T.peak_frequency, 'BinWidth', 0.5)
title('peak frequency')

subplot(1,2,2)
histogram(T.peak_SNR, 'BinWidth', 2)
title('SNR at peak frequency')


%% Optimized filter values

figure

subplot(1,4,1)
histogram(T.optim_window, 'Normalization', 'probability', 'BinWidth', 25)
title(['window_length [ ' sprintf('%i ', prctile(T.optim_window, [25 75])) ']'], 'Interpreter', 'none')

subplot(1,4,2)
histogram(T.optim_filter_ord, 'Normalization', 'probability', 'BinWidth', 15)
title(['filter_order [ ' sprintf('%i ', prctile(T.optim_filter_ord, [25 75])) ']'], 'Interpreter', 'none')

subplot(1,4,3)
histogram(T.optim_edge, 'Normalization', 'probability', 'BinWidth', 10)
title(['edge [ ' sprintf('%i ', prctile(T.optim_edge, [25 75])) ']'], 'Interpreter', 'none')

subplot(1,4,4)
histogram(T.optim_ar_ord, 'Normalization', 'probability', 'BinWidth', 5)
title(['ar_order [ ' sprintf('%i ', prctile(T.optim_ar_ord, [25 75])) ']'], 'Interpreter', 'none')

fprintf('\nMedian optimized parameters:\n')
default_parameters

%% Performance measures

phase_error_default = ang_diff(T.epochs_truephase_mean, T.epochs_estphase_default)';
phase_error_optim = ang_diff(T.epochs_truephase_mean, T.epochs_estphase_optim)';

figure
subplot(1,2,1)
histogram(rad2deg(ang_mean(phase_error_default)), 'BinWidth', 2)
hold on
histogram(rad2deg(ang_mean(phase_error_optim)), 'BinWidth', 2)

subplot(1,2,2)
histogram(rad2deg(ang_var2dev(ang_var(phase_error_default))), 'BinWidth', 3)
hold on
histogram(rad2deg(ang_var2dev(ang_var(phase_error_optim))), 'BinWidth', 3)

figure
polarhistogram(phase_error_default(:), 'NumBins', 120)
hold on
polarhistogram(phase_error_optim(:), 'NumBins', 120)


%% Effect of artifact rejection methods

%TODO


%%  SNR and angular variance of different estimation methods

figure
mdl = fitlm(T, 'error_ang_var_zrenner2018 ~ peak_SNR');
plot(mdl), hold on
mdl = fitlm(T, 'error_ang_var_default ~ peak_SNR');
plot(mdl), hold on
mdl = fitlm(T, 'error_ang_var_optim ~ peak_SNR');
plot(mdl)
ylabel('angular variance')


%% SNR and angular deviation of truephase (median, each epoch has a variance)

% calculate median angular deviation of truephase
T_calc = rowfun(@median, T(:, strcmp(T.Properties.VariableNames, 'epochs_truephase_ang_var')), 'OutputVariableNames', 'median_truephase_ang_var');
T_calc = rowfun(@(x) rad2deg(ang_var2dev(x)), T_calc, 'OutputVariableNames', 'median_truephase_ang_dev_deg');
T_calc.peak_SNR = T.peak_SNR;

figure
mdl = fitlm(T_calc, 'median_truephase_ang_dev_deg ~ peak_SNR');
plot(mdl), hold on

title('SNR and angular deviation of equivalent true phase estimates')
xlabel('peak SNR (dB)')
ylabel('circular deviation (deg)')

set(gcf, 'Renderer','Painters') %export vectorized
set(gcf, 'PaperPositionMode','manual','PaperType','a5','PaperOrientation','landscape','PaperUnits','normalized','PaperPosition',[0 0 1 1])
set(gcf, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(gcf, 'figure_snr_truephase_angdev', '-dpdf', '-r0')

clear('T_calc', 'mdl');


%%  SNR and angular deviation of truephase (median, each epoch has a variance) with error angular deviation (across epochs)
%(in units of circular variance, ranging from 0 .. 1)

% calculate median angular deviation of truephase
T_calc = rowfun(@median, T(:, strcmp(T.Properties.VariableNames, 'epochs_truephase_ang_var')), 'OutputVariableNames', 'error_truephase');
T_calc.peak_SNR = T.peak_SNR;
T_calc.error_optim = T.error_ang_var_optim;
T_calc.error_default = T.error_ang_var_default;
T_calc.error_zrenner2018 = T.error_ang_var_zrenner2018;

figure
mdl = fitlm(T_calc, 'error_truephase ~ peak_SNR');
plot(mdl), hold on
mdl = fitlm(T_calc, 'error_zrenner2018 ~ peak_SNR');
plot(mdl), hold on
mdl = fitlm(T_calc, 'error_default ~ peak_SNR');
plot(mdl), hold on
mdl = fitlm(T_calc, 'error_optim ~ peak_SNR');
plot(mdl)
ylabel('circular variance')

clear('T_calc', 'mdl');


%% SNR and angular deviation of truephase by epoch amplitude

%TODO

%% estimated amplitude vs true amplitude error

%TODO

%% true amplitude and phase error

amp_boundaries = [0 0.25 0.50 0.75 1.00];

figure
hold on

for i = 2:numel(amp_boundaries)

    T_calc = table();
    T_calc.amp = T.epochs_trueamp_mean;
    T_calc = [T_calc rowfun(@(x) quantile(x, amp_boundaries(i-1:i)), T_calc, 'OutputVariableNames', 'amp_boundaries')];
    T_calc.estphase_error = ang_diff(T.epochs_truephase_mean, T.epochs_estphase_default); % phase error
    T_calc = [T_calc rowfun(@(x, y, z) {z(x >= y(1) & x <= y(2))}, T_calc, 'OutputVariableNames', 'estphase_error_subset')]; % create new column with only included epochs
    T_calc = [T_calc rowfun(@(x) rad2deg(ang_var2dev(ang_var(x{1}))), T_calc(:,end), 'OutputVariableNames', 'estphase_error_subset_ang_dev')]; % calculate angular deviation
    T_calc.peak_SNR = T.peak_SNR;

    mdl = fitlm(T_calc, 'estphase_error_subset_ang_dev ~ peak_SNR');
    plot(mdl)

end

title('true amplitude and estimation error')


%% estimated amplitude and estimation error for all three methods

amp_boundaries = [0 0.25 0.50 0.75 1.00];

figure

algorithm_string = {'zrenner2018', 'default', 'optim'}
for j = 1:3
subplot(1,3,j)
    
color_order = get(0, 'DefaultAxesColorOrder');

plot_colors = {color_order(4,:), color_order(3,:), color_order(2,:), color_order(1,:)};
plot_handles = {};

for i = 1:(numel(amp_boundaries)-1)

    T_calc = table();
    T_calc.epochs_estamp = T.epochs_estamp_default;
    T_calc = [T_calc rowfun(@(x) quantile(x, amp_boundaries(i:(i+1))), T_calc, 'OutputVariableNames', 'amp_boundaries')];
    T_calc.estphase_error = ang_diff(T.epochs_truephase_mean, T.(['epochs_estphase_' algorithm_string{j}])); % phase error
    T_calc = [T_calc rowfun(@(x, y, z) {z(x >= y(1) & x <= y(2))}, T_calc, 'OutputVariableNames', 'estphase_error_subset')]; % create new column with only included epochs
    T_calc = [T_calc rowfun(@(x) rad2deg(ang_var2dev(ang_var(x{1}))), T_calc(:,end), 'OutputVariableNames', 'estphase_error_subset_ang_dev')]; % calculate angular deviation
    T_calc.peak_SNR = T.peak_SNR;

    mdl = fitlm(T_calc, 'estphase_error_subset_ang_dev ~ peak_SNR');
    plot_handles{i} = plot(mdl, 'Color', plot_colors{i}); hold on;

end

legend([plot_handles{1}(1), plot_handles{2}(1), plot_handles{3}(1), plot_handles{4}(1)], {'Lowest', 'Low', 'Medium', 'High'}, 'Location', 'southwest')

ylim([0 120])
    
title([algorithm_string{j} ' method estimation error by estimated amplitude quartile and SNR'])
xlabel('SNR (dB)')
ylabel('Estimation error circular deviation (deg)')

end
 

%% range or stdev and estimation error

%TODO: THIS DOESNT WORK

boundaries = [0 0.25 0.5 0.75 1.00];

figure
hold on

for i = 2:numel(boundaries)

    T_calc = table();
    T_calc.epochs_param = T.epochs_param_range_default;
    T_calc = [T_calc rowfun(@(x) quantile(x, boundaries(i-1:i)), T_calc, 'OutputVariableNames', 'param_boundaries')];
    T_calc.estphase_error = ang_diff(T.epochs_truephase_mean, T.epochs_estphase_default); % phase error
    T_calc = [T_calc rowfun(@(x, y, z) {z(x >= y(1) & x <= y(2))}, T_calc, 'OutputVariableNames', 'estphase_error_subset')]; % create new column with only included epochs
    T_calc = [T_calc rowfun(@(x) rad2deg(ang_var2dev(ang_var(x{1}))), T_calc(:,end), 'OutputVariableNames', 'estphase_error_subset_ang_dev')]; % calculate angular deviation
    T_calc.peak_SNR = T.peak_SNR;

    mdl = fitlm(T_calc, 'estphase_error_subset_ang_dev ~ peak_SNR');
    
    subplot(numel(boundaries), 2, (i-1)*2+1)
    plot(mdl)
    
    subplot(numel(boundaries), 2, (i-1)*2+2)
    polarhistogram(horzcat(T_calc.estphase_error_subset{:}))
    title(num2str(boundaries(i-1:i)))
end



%%


figure
mdl = fitlm(T, 'peak_SNR ~ peak_frequency');
plot(mdl)

figure
mdl = fitlm(T, 'optim_filter_ord ~ optim_window');
plot(mdl)


%does range predict phase error?
%range vs. epochs_estphase_optim

%does range predict phase estimatability?
%range vs. truephase_angdev

%sensitivity to parameter changes ?


%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  SNR and angular deviation of truephase (median, each epoch has a variance) with error angular deviation (across epochs)
%(in units of angular deviation)

% ----- FIGURE -----

figure_filename = 'fig_compare_parametersets';

% calculate median angular deviation of truephase
T_calc = rowfun(@median, T(:, strcmp(T.Properties.VariableNames, 'epochs_truephase_ang_var')), 'OutputVariableNames', 'median_truephase_ang_var');
T_calc = rowfun(@(x) rad2deg(ang_var2dev(x)), T_calc, 'OutputVariableNames', 'error_truephase');
T_calc.peak_SNR = T.peak_SNR;
T_calc.error_optim = rad2deg(ang_var2dev(T.error_ang_var_optim));
T_calc.error_default = rad2deg(ang_var2dev(T.error_ang_var_default));
T_calc.error_default_individualpassband = rad2deg(ang_var2dev(T.error_ang_var_default_individualpassband));
T_calc.error_zrenner2018 = rad2deg(ang_var2dev(T.error_ang_var_zrenner2018));
T_calc.error_zrenner2018_individualpassband = rad2deg(ang_var2dev(T.error_ang_var_zrenner2018_individualpassband));

figure('Color', 'w')

subplot(2, 3, [1 2 4 5])
hAx = gca; hold on

mdl = fitlm(T_calc, 'error_zrenner2018 ~ peak_SNR');
p1 = plot(mdl);
[p1.Color] = deal(hAx.ColorOrder(3,:));
[p1(2).LineWidth] = deal(2);
[p1(3:end).LineStyle] = deal('none');

mdl = fitlm(T_calc, 'error_zrenner2018_individualpassband ~ peak_SNR');
p2 = plot(mdl);
[p2.Color] = deal(hAx.ColorOrder(3,:));
[p2(2).LineStyle] = deal(':');
[p2(2).LineWidth] = deal(2);
[p2([3:end]).LineStyle] = deal('none');
p2(1).Marker = 'none';

mdl = fitlm(T_calc, 'error_default ~ peak_SNR');
p3 = plot(mdl);
[p3.Color] = deal(hAx.ColorOrder(2,:));
[p3(2).LineWidth] = deal(2);
[p3(3:end).LineStyle] = deal('none');

mdl = fitlm(T_calc, 'error_default_individualpassband ~ peak_SNR');
p4 = plot(mdl);
[p4.Color] = deal(hAx.ColorOrder(2,:));
[p4(2).LineStyle] = deal(':');
[p4(2).LineWidth] = deal(2);
[p4(3:end).LineStyle] = deal('none');

mdl = fitlm(T_calc, 'error_optim ~ peak_SNR');
p4(1).Marker = 'none';
p5 = plot(mdl);
[p5.Color] = deal(hAx.ColorOrder(1,:));
[p5(2).LineStyle] = deal(':');
[p5(2).LineWidth] = deal(2);
[p5(3:end).LineStyle] = deal('none');

mdl = fitlm(T_calc, 'error_truephase ~ peak_SNR');
p6 = plot(mdl)
[p6.Color] = deal([0.5 0.5 0.5]);
[p6(2).LineWidth] = deal(2);
[p6(3:end).LineStyle] = deal('none');

t1 = title('error deviation of different estimation methods')
xlabel('SNR (dB)')
ylabel('circular deviation (deg)')
legend([p1(2), p2(2), p3(2), p4(2), p5(2), p6(2)], {'500ms Window', 'with Individual Passband', '719ms Window', 'with Individual Passband', 'Individually Optimized', 'True Phase Error'})

subplot(2, 3, 3), hold on
histogram(T_calc.error_zrenner2018-T_calc.error_default, 'FaceColor', [0.7 0.7 0.7], 'BinWidth', 2)
histogram(T_calc.error_zrenner2018_individualpassband-T_calc.error_default_individualpassband, 'FaceColor', [0.4 0.4 0.4], 'BinWidth', 2)
legend({'8-13 Hz', 'Peak ±1Hz'})
t2 = title('Zrenner2018 vs. 719ms Window')
xlabel('improvement (deg)')
xlim([-4 30])

h2 = subplot(2, 3, 6), hold on
histogram(T_calc.error_zrenner2018-T_calc.error_zrenner2018_individualpassband, 'FaceColor', p1(1).Color, 'BinWidth', 2)
histogram(T_calc.error_default-T_calc.error_default_individualpassband, 'FaceColor', p3(1).Color, 'BinWidth', 2)
legend({'500ms Window', '719ms Window'})
t3 = title('Zrenner2018 vs. 719ms Window with Individual Passband')
xlabel('improvement (deg)')
xlim([-4 30])
ylim([0 80])

t1.String = 'A'; t1.FontSize = 14; t1.Position(1) = 1; t1.Position(2) = t1.Position(2) + 1;
t2.String = 'B'; t2.FontSize = 14; t2.Position(1) = 0; t2.Position(2) = t2.Position(2) + 1;
t3.String = 'C'; t3.FontSize = 14; t3.Position(1) = 0; t3.Position(2) = t3.Position(2) + 1;

set(gcf, 'Renderer','Painters') %export vectorized
set(gcf, 'PaperUnits', 'centimeter', 'PaperSize', [25 15]) % set size
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'normalized', 'PaperPosition',[0 0 1 1]); % fill page
set(gcf, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(gcf, figure_filename, '-dpdf', '-r0')

clear('T_calc', 'mdl');

%% true amplitude and estimation limit

% ----- FIGURE -----

figure_filename = 'fig_amp_true';
figure_textoutput_file_id = fopen([figure_filename '.txt'],'w');

amp_boundaries = [0 0.25 0.50 0.75 1.00];

figure('Color', 'w')
ax1 = subplot('Position',  [0.025 0.1 0.3 .825])
t1 = title('Filter Magnitude Response')

ax2 = subplot('Position',  [0.35 0.1 0.3 .825])

% calculate median angular deviation of truephase
T_calc = rowfun(@median, T(:, strcmp(T.Properties.VariableNames, 'epochs_truephase_ang_var')), 'OutputVariableNames', 'median_truephase_ang_var');
T_calc = rowfun(@(x) rad2deg(ang_var2dev(x)), T_calc, 'OutputVariableNames', 'median_truephase_ang_dev_deg');
T_calc.peak_SNR = T.peak_SNR;

mdl = fitlm(T_calc, 'median_truephase_ang_dev_deg ~ peak_SNR');
p = plot(mdl), hold on
[p.Color] = deal([0.5 0.5 0.5]);
[p(2).LineWidth] = deal(2);
    
fprintf(figure_textoutput_file_id, '\nrsquared=%1.4f, p=%1.1e', mdl.Rsquared.Adjusted, mdl.coefTest);

legend off, box on

t2 = title('SNR and angular deviation of equivalent true phase estimates')

xlabel('SNR (dB)')
ylabel('circular deviation (deg)')

ax3 = subplot('Position',  [0.675 0.1 0.3 .825])
hold on

color_order = get(0, 'DefaultAxesColorOrder');

plot_colors = {color_order(4,:), color_order(3,:), color_order(2,:), color_order(1,:)};
plot_handles = {};

for i = 1:numel(amp_boundaries)-1

    T_calc = table();
    T_calc.amp = T.epochs_trueamp_mean;
    T_calc = [T_calc rowfun(@(x) quantile(x, amp_boundaries(i:(i+1))), T_calc, 'OutputVariableNames', 'amp_boundaries')];
    T_calc.epochs_truephase_angdev = rad2deg(ang_var2dev(T.epochs_truephase_ang_var));
    T_calc = [T_calc rowfun(@(x, y, z) {z(x >= y(1) & x <= y(2))}, T_calc, 'OutputVariableNames', 'epochs_truephase_angdev_subset')]; % create new column with only included epochs
    T_calc = [T_calc rowfun(@(x) median(x{1}), T_calc(:,end), 'OutputVariableNames', 'epochs_truephase_angdev_subset_median')];
    T_calc.peak_SNR = T.peak_SNR;

    mdl = fitlm(T_calc, 'epochs_truephase_angdev_subset_median ~ peak_SNR');
    plot_handles{i} = plot(mdl, 'Color', plot_colors{i});
    [plot_handles{i}.Color] = deal(plot_colors{i});
    [plot_handles{i}(2).LineWidth] = deal(2);

end

legend([plot_handles{1}(1), plot_handles{2}(1), plot_handles{3}(1), plot_handles{4}(1)], {'Lowest', 'Low', 'Medium', 'High'})
box on

t3 = title('True estimation error by true amplitude quartile and SNR')

ylabel('Circular deviation of estimate (deg)')
ylabel('')
xlabel('SNR (dB)')

linkaxes([ax2 ax3], 'xy')
drawnow

t1.String = 'A';
t1.Position(1) = 1; t1.Position(2) = t1.Position(2) + 1;
t2.String = 'B';
t2.Position(1) = 1; t2.Position(2) = t2.Position(2) + 1;
t3.String = 'C';
t3.Position(1) = 1; t3.Position(2) = t3.Position(2) + 1;

set(gcf, 'Renderer','Painters') %export vectorized
set(gcf, 'PaperUnits', 'centimeter', 'PaperSize', [25 12]) % set size
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'normalized', 'PaperPosition',[0 0 1 1]); % fill page
set(gcf, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(gcf, figure_filename, '-dpdf', '-r0')

fclose(figure_textoutput_file_id);
clear('figure_filename', 'figure_textoutput', 'amp_boundaries', 'T_calc', 'mdl', 't1', 't2');


%% Estimated amplitude and estimation error for default method only

% ----- FIGURE -----

figure_filename = 'fig_amp_estimated';
figure_textoutput_file_id = fopen([figure_filename '.txt'],'w');

amp_boundaries = [0 0.25 0.50 0.75 1.00];

figure('Color', 'w')
axH_main = subplot(4,5,[1 2 3 4  6 7 8 9  11 12 13 14  16 17 18 19])
axH_main.Position = [0.075 0.08 0.725 .90];
axH_hist = {};
axH_hist{1} = subplot(4,5, 5, polaraxes); axH_hist{1}.Position = [0.79 0.77 0.21 .21];
axH_hist{2} = subplot(4,5,10, polaraxes); axH_hist{2}.Position = [0.79 0.54 0.21 .21];
axH_hist{3} = subplot(4,5,15, polaraxes); axH_hist{3}.Position = [0.79 0.31 0.21 .21];
axH_hist{4} = subplot(4,5,20, polaraxes); axH_hist{4}.Position = [0.79 0.08 0.21 .21];

color_order = get(0, 'DefaultAxesColorOrder');

plot_colors = {color_order(4,:), color_order(3,:), color_order(2,:), color_order(1,:)};
plot_handles = {};

for i = 1:(numel(amp_boundaries)-1)
    fprintf(figure_textoutput_file_id, 'amp range [%2.1f %2.1f]', amp_boundaries(i), amp_boundaries(i+1));
    
    T_calc = table();
    T_calc.epochs_estamp = T.epochs_estamp_default;
    T_calc = [T_calc rowfun(@(x) quantile(x, amp_boundaries(i:(i+1))), T_calc, 'OutputVariableNames', 'amp_boundaries')];
    T_calc.estphase_error = ang_diff(T.epochs_truephase_mean, T.epochs_estphase_default); % phase error
    T_calc = [T_calc rowfun(@(x, y, z) {z(x >= y(1) & x <= y(2))}, T_calc, 'OutputVariableNames', 'estphase_error_subset')]; % create new column with only included epochs
    T_calc = [T_calc rowfun(@(x) rad2deg(ang_var2dev(ang_var(x{1}))), T_calc(:,end), 'OutputVariableNames', 'estphase_error_subset_ang_dev')]; % calculate angular deviation
    T_calc.peak_SNR = T.peak_SNR;

    mdl = fitlm(T_calc, 'estphase_error_subset_ang_dev ~ peak_SNR');
    axes(axH_main)
    plot_handles{i} = plot(mdl); hold on;
    [plot_handles{i}.Color] = deal(plot_colors{i});
    [plot_handles{i}(2).LineWidth] = deal(2);
    
    fprintf(figure_textoutput_file_id, ' - rsquared=%1.4f, p=%1.1e', mdl.Rsquared.Adjusted, mdl.coefTest);
    
    axes(axH_hist{i})
    polarhistogram(horzcat(T_calc.estphase_error_subset{:}), 'Normalization', 'probability', 'BinWidth', pi/18, 'FaceColor', plot_colors{i})
    set(gca, 'RTick', [])
    set(gca, 'ThetaTickLabel', [])
    set(gca, 'GridColor', [0 0 0], 'GridAlpha', 1)

    fprintf(figure_textoutput_file_id, ' - error circ deviation %1.1f', rad2deg(ang_var2dev(ang_var(horzcat(T_calc.estphase_error_subset{:})))));
    fprintf(figure_textoutput_file_id, '\n');
end

axes(axH_main)
legend(axH_main, [plot_handles{1}(1), plot_handles{2}(1), plot_handles{3}(1), plot_handles{4}(1)], {'Lowest', 'Low', 'Medium', 'High'}, 'Location', 'northeast')

title('Estimation error by estimated amplitude quartile and SNR')
title('')
xlabel('SNR (dB)')
ylabel('circular deviation (deg)')

set(gcf, 'Renderer','Painters') %export vectorized
set(gcf, 'PaperPositionMode','manual','PaperType','a5','PaperOrientation','landscape','PaperUnits','normalized','PaperPosition',[0 0 1 1])
set(gcf, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(gcf, figure_filename, '-dpdf', '-r0')

clear('amp_boundaries', 'T_calc', 'mdl');


%% Visualize Filter Objects

% ----- FIGURE -----

filter_objects = {};
filter_labels = {};
filter_colorindex = [];
peak_frequency = 11;
fs = 1000;
%FIR
for ord = [2 3 4 5] % multiple of the period, notice that the filters are applied forward and backward    
    % windowed sinc
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'CutoffFrequency1', peak_frequency-1, 'CutoffFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'window')};
    filter_labels = {filter_labels{:} ['window ' num2str(ord) 'x']};
    filter_colorindex = [filter_colorindex 1];
end
for ord = [3 4 5] % multiple of the period, notice that the filters are applied forward and backward    
    % least squares (equiripple is similar)
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'StopbandFrequency1', peak_frequency-4, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+4, 'SampleRate', fs, 'DesignMethod', 'ls')};
    filter_labels = {filter_labels{:} ['ls ' num2str(ord) 'x']};
    filter_colorindex = [filter_colorindex 2];
end
%IIR
for ord = [4 8 12]
    % butterworth
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'HalfPowerFrequency1', peak_frequency-1, 'HalfPowerFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'butter')};
    filter_labels = {filter_labels{:} ['butter ' num2str(ord)]};
    filter_colorindex = [filter_colorindex 3];    
end
for ord = [4 6 8]
    % chebychev I
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'PassbandFrequency1', peak_frequency-1.5, 'PassbandFrequency2', peak_frequency+1.5, 'PassbandRipple', 0.5, 'SampleRate', fs, 'DesignMethod', 'cheby1')};
    filter_labels = {filter_labels{:} ['cheby1 ' num2str(ord)]};
    filter_colorindex = [filter_colorindex 4];        
end
for attenuation = [10 20]
    % elliptic
    filter_objects = {filter_objects{:} designfilt('bandpassiir', 'StopbandFrequency1', peak_frequency-3, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+3, 'StopbandAttenuation1', attenuation, 'PassbandRipple', 0.5, 'StopbandAttenuation2', attenuation, 'SampleRate', fs, 'DesignMethod', 'ellip', 'MatchExactly', 'passband')};
    filter_labels = {filter_labels{:} ['ellip ' num2str(2*attenuation) 'dB']};
    filter_colorindex = [filter_colorindex 5];
end

fprintf('\n%i filter objects created\n', numel(filter_objects));

h = fvtool(filter_objects{:});
title('')

h.zoom('x', [0 35]);
h.zoom('y', [-40 2.5]); % labels will be doubled later
legend(h, filter_labels{:}, 'Location', 'northeast');

%Note: we need to square the amplitude response to show the effect of
%forward and backward filtering, i.e. double the logarithm dB response
h.CurrentAxes.YTick = [-50:5:0];
h.CurrentAxes.YTickLabel = string(h.CurrentAxes.YTick * 2);

for i = 1:numel(filter_colorindex)
    h.CurrentAxes.Children(end-i+1).Color = h.CurrentAxes.ColorOrder(filter_colorindex(i),:);
end

ax = h.Children(end);

f = figure;
ax.Parent = f;



set(h,'Renderer','Painters') %export vectorized
set(h,'PaperPositionMode','manual','PaperType','a5','PaperOrientation','landscape','PaperUnits','normalized','PaperPosition',[0 0 1 1])
set(h, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(h,'figure_filter_visualization', '-dpdf', '-r0')


fprintf('\nAngular deviation sanity check:')
[ang_m, ang_v, amp_m, amp_v] = phastimate_truephase(cos((1:(2*fs))./fs*2*pi*11)', filter_objects);
fprintf('\nSignal only:    %5.1f deg', rad2deg(sqrt(2*ang_v)))
[ang_m, ang_v, amp_m, amp_v] = phastimate_truephase(cos((1:(2*fs))./fs*2*pi*11)' + randn(2000,1)*10, filter_objects);
fprintf('\nSignal + Noise: %5.1f deg', rad2deg(sqrt(2*ang_v)))
[ang_m, ang_v, amp_m, amp_v] = phastimate_truephase(randn(2000,1), filter_objects);
fprintf('\nNoise only:     %5.1f deg', rad2deg(sqrt(2*ang_v)))
fprintf('\n')
