function [D_out, cfg_out, diff_sig_mask] = hpf_and_instantaneous_artifact_mask(D, cfg)
% Inputs: D, spike_dur, iqr_thr
% Outputs: D
% Description: Identifies instantaneous artifacts and masks them. Then performs HPF

vardefault('cfg',struct);
field_default('cfg', 'spike_dur', 0.1) % estimated duration of a spike artifact in seconds, from spike onset to peak
field_default('cfg', 'iqr_thr', 3) % threshold to identify outliers (e.g. outlier > 75th percentile + iqr_thr*interquartile range)
field_default('cfg', 'f_c', 2) % % Cutoff frequency for high-pass filter

D_out = D; 

og_sig = cat(2, D_out.trial{:});                             % og_sig contains an electrodes*timepoints matrix of "original values"

% IDENTIFY POTENTIAL OUTLIERS (as samples where diff_sig > Q3 + n_thr * IQR, or diff_sig < Q1 - n_thr * IQR)
diff_sig = diff(og_sig,1,2);                                % diff_sig contains an electrodes*(timepoints-1) matrix of differences (temporal derivative)
m = 2*round(cfg.spike_dur*D_out.fsample/2) + 1; % force m to be odd
diff_sig_smoothed = convn(diff_sig(:,max(1,min(size(diff_sig,2),1-(m-1)/2:size(diff_sig,2)+(m-1)/2))), hanning(m)', 'valid');

% [iqr_diff,qart_diff] = iqr(diff_sig_smoothed,2);   % iqr_diff: interquartile range of differences (IQR); qart_diff: first and third quartiles (Q1 & Q3) (REQUIRES >=R2024a)
iqr_diff = iqr(diff_sig_smoothed,2);
qart_diff = prctile(diff_sig_smoothed, [25; 75], 2);

% RECONSTRUCTED SIGNAL
diff_sig_mask = diff_sig_smoothed > qart_diff(:,2)+cfg.iqr_thr*iqr_diff | diff_sig_smoothed < qart_diff(:,1)-cfg.iqr_thr*iqr_diff; % crops derivatives beyond minimum/maximum values
diff_sig(diff_sig_mask) = 0;

og_sig = cumsum([zeros(size(diff_sig,1),1), diff_sig],2);   % reconstructs original signal by cumulative sum (temporal integral)

% HPF
k = 1/(1 + 2*pi*cfg.f_c/D_out.fsample); % first order IIR
for n=2:size(diff_sig,2)
    og_sig(:,n) = k*og_sig(:,n-1) + k*diff_sig(:,n-1); % HPF
end

D_out.trial = {og_sig};
cfg_out = cfg; 