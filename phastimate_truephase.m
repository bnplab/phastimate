function [ang_m, ang_v, amp_m, amp_v] = phastimate_truephase(data, filter_objects)
%PHASTIMATE_TRUEPHASE Estimate phase at center of epoch using multiple equivalent methods
%   [ang_m, ang_v, amp_m, amp_v] = phastimate_truephase(data, filter_objects)
%
%   Input:
%     data is a time x epoch matrix
%     filter_objects is a cell array of digitalFilter
%
%   Output:
%     ang_m  circular mean of phase as estimated with different filters
%     ang_v  circular variance of phase
%     amp_m  mean of amplitude
%     amp_v  variance of amplitude
%
%   Example:
%     Fs = 200;
%     data = [sin([1:Fs]/Fs*5*2*pi)' sin([1:Fs]/Fs*5.5*2*pi)' sin([1:Fs]/Fs*6*2*pi)'] + randn(Fs,3);
%     D1 = designfilt('bandpassfir', 'FilterOrder', round(Fs/5), 'CutoffFrequency1', 4, 'CutoffFrequency2', 7, 'SampleRate', Fs);
%     D2 = designfilt('bandpassiir', 'FilterOrder', 6, 'HalfPowerFrequency1', 3, 'HalfPowerFrequency2', 8, 'SampleRate', Fs);
%
%     [ang_m, ang_v, amp_m, amp_v] = phastimate_truephase(data, {D1 D2})

assert(iscell(filter_objects), 'filter_objects must be a cell array of digitalFilter objects')
cellfun(@(x) assert(isa(x, 'digitalFilter'), 'filter_objects cell array members must be of type digitalFilter'), filter_objects)

% demean the data
data = detrend(data,'constant');

% filter the data through each of the filter objects and apply Hilbert transform
data = cellfun(@(x) hilbert(filtfilt(x, data)), filter_objects, 'UniformOutput', false);

% take the value at the window center of each epoch
analytic = cell2mat(cellfun(@(x) x(ceil(end/2),:), data, 'UniformOutput', false)');

%visualize for debugging
%figure, plot(real(analytic), imag(analytic))

% phase
% calculate the mean normalized phasor
mean_normalized = mean(analytic ./ abs(analytic));

ang_m = angle(mean_normalized); % circular mean
ang_v = 1-abs(mean_normalized); % circular variance

% amplitude
amp_m = mean(abs(analytic));
amp_v = var(abs(analytic));

end

