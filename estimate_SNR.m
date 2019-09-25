%% estimate SNR
% cut the signal into 4s snippets for artefact rejection

function [peak_frequency, peak_SNR] = estimate_SNR(epochs, fs, frequency_interval, ax)


[nr_samples, nr_trials] = size(epochs);
segment_length = 2*fs;

[power,frequency] = pwelch(epochs(:),hann(segment_length),0, nr_samples,fs);

bins_selected = find(frequency>=2,1):find(frequency>=45,1);
frequency = frequency(bins_selected,:)';
power = power(bins_selected,:);

log_freq = log10(frequency);
log_power = 10*log10(power);

% define frequence bands
idx1 = find(frequency>=0.5 & frequency <= 7);
idx2 = find(frequency>=35 & frequency <= 65);
idx3 = find(frequency>=frequency_interval(1) & frequency <= frequency_interval(2));
idx4 = find(frequency>=.5 & frequency<65);

idx_freq = [idx1 idx2];
log_freq_roi = [ones(1,numel(idx_freq)); log_freq(idx_freq)];
log_power_roi = [log_power(idx1)' log_power(idx2)'];


% fit linear slope in log-log space

fit = log_freq_roi'\log_power_roi';
slope = fit(2);
intercept = fit(1);
fit_1f = slope*log10(frequency') + intercept;
power_corrected = log_power - fit_1f;

% find alpha peak
[pks, locs] = findpeaks(power_corrected, 'MinPeakProminence',0.5);
locs_alpha = intersect(locs, idx3);
peak_SNR = power_corrected(locs_alpha);
noise_level = mean(power_corrected(idx1));
spect = frequency(idx4);

% select the maximum alpha_peak
[~, idx_max] = max(peak_SNR);
peak_frequency = spect(locs_alpha(idx_max));
peak_SNR = peak_SNR(idx_max);

% plot
if ~isempty(ax) == 1
    
    semilogx(ax, frequency, log_power)
    hold on;
    semilogx(ax, frequency(idx_freq), log_power_roi, '.')
    semilogx(ax, frequency, fit_1f)
    semilogx(ax, frequency, power_corrected)
    
    semilogx(ax, peak_frequency, peak_SNR, 'rs')
    semilogx(ax, frequency(idx1), repmat(noise_level, [numel(idx1), 1]),'g')
    
    xlim([2 50])
    
    fprintf('peak at freq = %.2f Hz\n', peak_frequency)
    fprintf('    with SNR = %.2f db\n', peak_SNR)
    
    %title(sprintf('peak = %.3f, crit=%i', peak_SNR, peak_SNR > 10))
    if isempty(peak_frequency)
        title('no peak found')
    else
        title(sprintf('%.1f dB SNR at %.1f Hz', peak_SNR, peak_frequency))
    end
    legend( 'original spectrum', 'points used for fitting', '1/f fit', ...
            'adjusted spectrum')
    xlabel('frequency (Hz)')
    ylabel('power (dB)')
end


end
