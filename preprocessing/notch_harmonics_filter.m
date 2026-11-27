%%% run notch filtering to remove line noise - 60hz and harmonics
% bml repo no longer required on path

function [D_filt, cfg_out] = notch_harmonics_filter(D_in, cfg)

vardefault('cfg', struct);

% Determine sample rate (assume 1000 Hz if not specified on D_in)
if isfield(D_in, 'fsample') && ~isempty(D_in.fsample)
    Fs = D_in.fsample;
else
    Fs = 1000;
end
nyq = Fs / 2;

% modifiable parameters (clip default harmonics to stay below Nyquist)
default_harmonics = 60:60:480;
field_default('cfg', 'dftfreq',            default_harmonics(default_harmonics < (nyq - 2)));
field_default('cfg', 'notch_bandwidth_hz', 2);

% Also clip user-supplied dftfreq if any exceed Nyquist
cfg.dftfreq = cfg.dftfreq(cfg.dftfreq < (nyq - 2));

% fix the following parameters
cfg.dftfilter = 'yes';

% % % % % % % % % % % % % % 2026-7-1 note: this section originally was dealing with problem of not always having complete 60hz cycles
% % % % % % % % % % % % % %%%% RD not sure if even necessary anymore
% % % % % % % % % % % % % D_annot = bml_raw2annot(D_in);
% % % % % % % % % % % % % D_annot.nSamples2 = round(floor(D_annot.nSamples .* cfg.dftfreq(1)./D_annot.Fs) .* D_annot.Fs./cfg.dftfreq(1));
% % % % % % % % % % % % % D_annot.nSamples2 = round(floor(D_annot.nSamples2 .* cfg.dftfreq(1)./D_annot.Fs) .* D_annot.Fs./cfg.dftfreq(1));
% % % % % % % % % % % % % D_annot.nSamples2 = D_annot.nSamples2(:,1);
% % % % % % % % % % % % % D_annot.ends = D_annot.starts + D_annot.nSamples2 ./ D_annot.Fs;
% % % % % % % % % % % % % 
% % % % % % % % % % % % % cfg1=[];
% % % % % % % % % % % % % cfg1.epoch=D_annot;
% % % % % % % % % % % % % D1= bml_redefinetrial(cfg1,D_in);
% % % % % % % % % % % % % 
% % % % % % % % % % % % % mask zeroes as nans
% % % % % % % % % % % % cfg1=[];
% % % % % % % % % % % % cfg1.remask_nan = true;
% % % % % % % % % % % % cfg1.value = 0;
% % % % % % % % % % % % D1 = bml_mask(cfg1, D1);
% % % % % % % % % % % % F = cfg.dftfreq;
% % % % % % % % % % % % N = D_annot.nSamples2;
% % % % % % % % % % % % Fs = D_annot.Fs;

%%%%%%%%%% 2026/11/27 AM replaced above section with new section below to remove BML toolbox dependencies %%%%%%%
% Replicate bml_raw2annot + bml_redefinetrial + bml_mask without bml toolbox
D1 = D_in;
nTrials = numel(D_in.trial);
nSamples_per_trl = cellfun(@(x) size(x, 2), D_in.trial);

nSamples2_vec = round(floor(nSamples_per_trl .* cfg.dftfreq(1) ./ Fs) .* Fs ./ cfg.dftfreq(1));
nSamples2_vec = round(floor(nSamples2_vec .* cfg.dftfreq(1) ./ Fs) .* Fs ./ cfg.dftfreq(1));

for iTrl = 1:nTrials
    nSmp2 = nSamples2_vec(iTrl);
    trl_dat = D_in.trial{iTrl}(:, 1:nSmp2);
    % Replicate bml_mask (cfg1.remask_nan = true; cfg1.value = 0)
    trl_dat(trl_dat == 0) = NaN;
    D1.trial{iTrl} = trl_dat;
    D1.time{iTrl}  = D_in.time{iTrl}(1:nSmp2);
end

if isfield(D1, 'sampleinfo') && size(D1.sampleinfo, 1) == nTrials
    D1.sampleinfo(:, 2) = D1.sampleinfo(:, 1) + nSamples2_vec(:) - 1;
end

F = cfg.dftfreq;
N = nSamples2_vec(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newF = F + reshape(find(abs(1/N*Fs*(0:N-1) - F(1)) < cfg.notch_bandwidth_hz / 2) - 1, [], 1)/N*Fs - F(1);
newF = reshape(newF, 1, []);
cfg.dftfreq = newF;

% intended to use Rohan's modified version of fieldtrip functions; this will not work right if normal FT repo is higher on your path
D_filt = ft_preprocessing(cfg, D1);

cfg_out = cfg;