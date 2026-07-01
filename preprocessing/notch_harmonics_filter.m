%%% run notch filtering to remove line noise - 60hz and harmonics
% bml repo must be on path

function [D_filt, cfg_out] = notch_harmonics_filter(D_in, cfg)

vardefault('cfg',struct);

% modifiable parameters
field_default('cfg','dftfreq',          [60 120 180 240 300 360 420 480]);
field_default('cfg','dftbandwidth',     [1   1   1   1   1   1   1   1]);
field_default('cfg','dftneighbourwidth',[2   2   2   2   2   2   2   2]);

% fix the following parameters
cfg.dftfilter='yes';

D_annot = bml_raw2annot(D_in);
D_annot.nSamples2 = round(floor(D_annot.nSamples .* cfg.dftfreq(1)./D_annot.Fs) .* D_annot.Fs./cfg.dftfreq(1));
D_annot.nSamples2 = round(floor(D_annot.nSamples2 .* cfg.dftfreq(1)./D_annot.Fs) .* D_annot.Fs./cfg.dftfreq(1));
D_annot.nSamples2 = D_annot.nSamples2(:,1);


%% AM question  - does this still need to be here? should it be a conditional? 
% D_annot.nSamples2 = 938640; % hardcoded edit for DM1037


D_annot.ends = D_annot.starts + D_annot.nSamples2 ./ D_annot.Fs;

cfg1=[];
cfg1.epoch=D_annot;
D1= bml_redefinetrial(cfg1,D_in);

%% AM note: skipping selectdata step because this would be awkward to use when all chans are already plotted in gui, including audio

% mask zeroes as nans
cfg1=[];
cfg1.remask_nan = true;
cfg1.value = 0;
D1 = bml_mask(cfg1, D1);


%% am note: in clean_mask_hpf_notch_filter, this is where we'd run hpf and cleaning before notch filtering, but in this case the data is already cleaned+hpf'ed...
... is it ok that we're doing masking and then jumping straight to notch filtering? 


F = cfg.dftfreq;
N = D_annot.nSamples2;
Fs = D_annot.Fs;
BW = 2; % Bandwith of notch filter in Hz

newF = F + reshape(find(abs(1/N*Fs*(0:N-1) - F(1))<BW/2)-1,[],1)/N*Fs - F(1);
newF = reshape(newF,1,[]);
cfg.dftfreq = newF;

% intended to use Rohan's modified version of fieldtrip functions
D_filt = ft_preprocessing(cfg,D1);

cfg_out = cfg; 