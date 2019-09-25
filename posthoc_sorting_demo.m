%% Post-hoc sorting demo

% switch to current directory and add relative path to phastimate toolbox
cd(fileparts(getfield(matlab.desktop.editor.getActive, 'Filename')))
addpath('../phastimate_code/')

% ----- FIGURE -----

mepdata = load(fullfile('..', 'data', 'mepdata3.mat')); 

D = designfilt('bandpassfir', 'FilterOrder', 190, 'CutoffFrequency1', 8, 'CutoffFrequency2', 13, 'SampleRate', 1000, 'DesignMethod', 'window');
estphase = phastimate(double(mepdata.prestimulus_eeg(750:1499,:)), D, 64, 30, 128);

%wrap to -1.5*pi .. 0.5*pi
estphase(estphase > 0.5*pi) = estphase(estphase > 0.5*pi)-(2*pi);

y = log(mepdata.mep_amplitude(:,1));

mdl = fitlm([cos(estphase)' sin(estphase)'], y)

coeffs = mdl.Coefficients.Estimate;

%(intercept) + x1*cos(x) + x2*sin(x)
bestfit = @(x) coeffs(1) + coeffs(2)*cos(x) + coeffs(3)*sin(x);
phasefit = @(x) coeffs(1) + cos(x);
%phasefit = @(x) coeffs(1) + abs(cos(x*0.5))*4 - 2;

%% Power Analysis using bootstrapping
rng('shuffle')

nsample_range = 20:5:200;
nboot = 2000;
pvalue_contig = nan(numel(nsample_range), nboot);
pvalue_indep = nan(numel(nsample_range), nboot);

%TODO: use randsample with replacement?

for i = 1:numel(nsample_range)
    nsample = nsample_range(i)
    for j = 1:nboot
        %contiguous samples
        %indices = randsample(length(estphase)-nsample+1, 1) + (0:(nsample-1));
        %mdl = fitlm([cos(estphase(indices))' sin(estphase(indices))'], y(indices));
        %pvalue_contig(i, j) = mdl.coefTest;
        
        %independent samples
        indices = randsample(length(estphase), nsample);
        mdl = fitlm([cos(estphase(indices))' sin(estphase(indices))'], y(indices));
        pvalue_indep(i, j) = mdl.coefTest;
        
    end
end


%%
figure
mVticks = [0.05 0.1 0.2 0.5 1 2 5];
histogram(y, 'FaceColor', [0.5 0.5 0.5])
set(gca, 'XTick', log(mVticks*1e3), 'XTickLabel', string(num2cell(mVticks)))
xlabel('mV')

%%
figure('Color', 'white')

ax = subplot('Position', [0.05 0.2 0.6 0.75])

scatter(estphase, y, 'SizeData', 150, 'MarkerFaceColor', [0.6 0.6 0.7], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0), hold on
fplot(bestfit, 'r', [-1.5*pi 0.5*pi], 'LineWidth', 2.5)
fplot(phasefit, [-1.5*pi 0.5*pi], 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5, 'LineStyle', ':')
%TODO: plot CI (bootstrapping residuals?)
xlim([-pi pi])
ax.YTick = log(mVticks*1e3); ax.YTickLabel = string(num2cell(mVticks));
ylabel('amplitude (mV)')
title(sprintf('sine F-statistic vs. constant p-value = %1.3g', mdl.coefTest))
hold on
h = boxplot(y, discretize(estphase, (-1.5*pi):(pi/5):(0.5*pi)), 'positions', (-1.5*pi+pi/10):(pi/5):(0.5*pi-pi/10), 'width', pi/5, 'Whisker', 0, 'Symbol', '');
set(h,{'LineWidth'},{1.25})
ax.XTick = (-1.5*pi+pi/10):(pi/5):(0.5*pi-pi/10);
ax.XTickLabel = cellstr(num2str(rad2deg(ax.XTick'), '%i°'));
ax.XTickLabel(ax.XTick == -pi) = {'trough'};
ax.XTickLabel(ax.XTick == 0) = {'peak'};
ax.FontSize = 12; %ax.FontWeight = 'bold';
axis tight
xlabel('pre-stimulus phase')
t1 = title('A')

for i = 1:10
    subplot('Position', [0.05+0.1*0.6*(i-1)+0.015 0.025 0.08*0.6 0.08])
    plot(mean(mepdata.prestimulus_eeg(end-150:end,discretize(estphase, (-1.5*pi):(pi/5):(0.5*pi)) == i), 2), 'Color', 'b', 'LineWidth', 1.5)
    box off
    xlim([0 170]), ylim([-6.5 5])
    set(gca, 'Color', [1 0.9 0.9], 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none')
    rectangle('Position', [149 -6.5 21 11.5], 'LineStyle', 'none', 'FaceColor', [.5 .5 .5])
end

%
%figure('Color', 'white')
subplot('Position', [0.725 0.2 0.25 0.75])
%plot(nsample_range, mean(pvalue_contig < 0.05, 2), 'LineWidth', 2)
hold on
plot(nsample_range, mean(pvalue_indep < 0.05, 2), 'Color', 'b', 'LineWidth', 2)
plot([nsample_range(1) nsample_range(end)], [0.8 0.8], 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', ':')
box on
axis tight
ax = gca;
ax.FontSize = 12; %ax.FontWeight = 'bold';
xlabel('trials')
ylabel('power')
t2 = title('B')

t1.FontSize = 14; t1.Position(1) = -5.1;
t2.FontSize = 14; t2.Position(1) = 1

%%

set(gcf, 'Renderer','Painters') %export vectorized
set(gcf, 'PaperUnits', 'centimeter', 'PaperSize', [25 15]) % set size
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'normalized', 'PaperPosition',[0 0 1 1]); % fill page
set(gcf, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard
print(gcf, 'fig_post_hoc_sorting', '-dpdf', '-r0')