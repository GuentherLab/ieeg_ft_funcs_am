%%%% created trialwise and trial-averaged spectrograms for individual channels

op.sub = 'DBS3012';
op.art_crit = 'F'; 
op.tf_rate = 20; %Hz Sampling rate of time frequency plot
op.tf_foi = round(10.^(0.30:0.05:2.4),2,'signif');

set(0,'DefaultFigureWindowStyle','docked')

set_project_specific_variables()

ft_to_load = [PATH_FIELDTRIP, filesep, op.sub, '_ft_raw_filt_trial.mat'];

%% load data
load(ft_to_load)


% dbs and macro only
cfg = [];
cfg.channel = {'dbs*','macro*'};
D_sel = ft_selectdata(cfg,D);

%% need to change D_sel.time so that times are relative to a timepoint within each trial, like speech onset
%% then pick a window around this timepoint as input to ft_freqanalysis via cfg.toi


%% make spectrogram 
cfg=[];
cfg.foi = op.tf_foi;
cfg.dt = 1/op.tf_rate;
cfg.trials = [1 384]; 
cfg.output       = 'pow';
cfg.method       = 'wavelet';
cfg.width  = 7; % number of wavelets
cfg.feedback     = 'no';
cfg.keeptrials   = 'yes';
cfg.keeptapers   = 'no';
cfg.pad          = 'nextpow2';

%%%%%%%%%%%%%%% calcuations to get toi 
dt        = bml_getopt(cfg,'dt',0.020); %width of pixels in returned matrix
epoch     = bml_raw2annot(D_sel);
% cfg.toi = mean([epoch.starts, epoch.ends],2); 
cfg.toi = D_sel.time{386};


% % % tic;             tf = bml_freqanalysis_power_wavelet(cfg, D_sel);          toc
tic;             tf = ft_freqanalysis(cfg, D_sel);          toc

% plot sample spectrogram
chind = 1; 
rpt=2;
chdat = squeeze(tf.powspctrm(:,chind,:,:));
spgrm = squeeze(chdat(rpt,:,:)); 
imagesc(flipud(spgrm)) % flip to put low freqs on the bottom
hax = gca;

% make y axis labels
nyticks = 15; % number of yticks to use
nfreq = length(tf.freq);
yfreq_inds = round(linspace(1,nfreq,nyticks)); 
ytick_freqs = tf.freq(yfreq_inds); 
hax.YTick = yfreq_inds';
hax.YTickLabel = num2str(round(flipud(ytick_freqs'),2,'signif')); 



colorbar