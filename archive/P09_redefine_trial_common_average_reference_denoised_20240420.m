%% currently only working with Triplet

%%% AM; modified 2022/8/16 to work with denoised trialwise data

function P09_redefine_trial_common_average_reference_denoised(SUBJECT)

% % % % % Loading packages
ft_defaults
bml_defaults
format long

% % % % % Loading parameters
% SUBJECT='DBS3012';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_SAVE_PREPROCESSED = [PATH_SUBJECT '/Preprocessed Data/FieldTrip'];
cd(PATH_SYNC)

% padding time to include before/after trial onset/offset
prebuffer_sec = 1.5; 
postbuffer_sec = 2; 

% % % Load annotation tables

electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
artifact = bml_annot_read(['annot/' SUBJECT '_artifact_criteria_E_denoised.txt']); % use multi-frequency-averaged high gamma
empty_electrode = bml_annot_read(['annot/' SUBJECT '_empty_electrode.txt']);
cue_presentation = bml_annot_read(['annot/' SUBJECT '_cue_presentation.txt']);

%%
% % % Load FieldTrip raw data

load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_raw_filt_trial_denoised.mat']);

% % % Adjusting length of sessions for notch filter to work
% % % % % % D_annot = bml_raw2annot(D);
% % % % % % Fl=[60 120 180 240 300 360 420 480];
% % % % % % D_annot.nSamples2 = round(floor(D_annot.nSamples .* Fl(1)./D_annot.Fs) .* D_annot.Fs./Fl(1));
% % % % % % D_annot.nSamples2 = round(floor(D_annot.nSamples2 .* Fl(1)./D_annot.Fs) .* D_annot.Fs./Fl(1));
% % % % % % D_annot.nSamples2 = D_annot.nSamples2(:,1);
% % % % % % D_annot.ends = D_annot.starts + D_annot.nSamples2 ./ D_annot.Fs;
% % % % % % 
% % % % % % cfg=[];
% % % % % % cfg.epoch=D_annot;
% % % % % % D1= bml_redefinetrial(cfg,D);


% % % Selecting electrodes and remasking with zeros instead of NaNs

% % % Some of FieldTrip's functions don't work with NaNs, so we are going to temporary replace NaNs with zeros to avoid issues. 

cfg=[];
cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'};
% % %         cfg.trials = logical([1 1 1 1 0]);
D_sel = ft_selectdata(cfg,D);

% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % cfg=[];
% % % % % % % % % % % % % % % % % % cfg.remask_nan = true;
% % % % % % % % % % % % % % % % % % cfg.value = 0;
% % % % % % % % % % % % % % % % % % D_sel = bml_mask(cfg, D_sel);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % Applying high pass filter and line noise removal filter
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % cfg=[];
% % % % % % % % % % % % % % % % % % cfg.hpfilter='yes';
% % % % % % % % % % % % % % % % % % cfg.hpfreq=1;
% % % % % % % % % % % % % % % % % % cfg.hpfilttype='but';
% % % % % % % % % % % % % % % % % % cfg.hpfiltord=5;
% % % % % % % % % % % % % % % % % % cfg.hpfiltdir='twopass';
% % % % % % % % % % % % % % % % % % cfg.dftfilter='yes';
% % % % % % % % % % % % % % % % % % %%%% this interpolation option is causing an error with 3012, artifact criterion E
% % % % % % % % % % % % % % % % % % % cfg.dftreplace='neighbour';  %using spectrum interpolation method Mewett et al 2004
% % % % % % % % % % % % % % % % % % cfg.dftfreq           = [60 120 180 240 300 360 420 480];
% % % % % % % % % % % % % % % % % % cfg.dftbandwidth      = [ 1   1   1   1   1   1   1   1];
% % % % % % % % % % % % % % % % % % cfg.dftneighbourwidth = [ 2   2   2   2   2   2   2   2];
% % % % % % % % % % % % % % % % % % D_sel_filt = ft_preprocessing(cfg,D_sel);

%% Redefining trials

% % % If the cue_presentation table has other column names, modify the following lines accordingly. 
% % % Use generous buffers around the cue and speech production, even if this means that consecutive trials overlap. 

epoch = cue_presentation(:,{'stim1_onset','ends','session_id','trial_id'});
epoch.starts = epoch.stim1_onset - prebuffer_sec;
epoch.ends = epoch.ends + postbuffer_sec;
epoch = bml_annot_table(epoch);

cfg = [];
cfg.epoch = epoch;
D_sel_trial = bml_redefinetrial(cfg,D_sel); % D was already in trial format, but we are redefining trial timing


% % % Masking artifacts and empty electrodes with NaNs for re-referencing

%combining artifacts with empty_electrode table
cfg=[];
cfg.groupby='label';
artifact_empty = bml_annot_union(cfg, artifact, empty_electrode);

%masking artifacts and empty_electrodes with NaNs
cfg=[];
cfg.annot=artifact_empty;
cfg.label_colname = 'label';
cfg.complete_trial = true; %masks entire trials
cfg.value=NaN;
D_sel_trial_mask = bml_mask(cfg,D_sel_trial);

% % % Common trimmed average reference per connector groups

el_ecog = electrode(electrode.type=="ecog",:);

cfg=[];
cfg.label = el_ecog.electrode;
cfg.group = el_ecog.connector;
cfg.method = 'CTAR'; %using trimmed average referencing
cfg.percent = 50; %percentage of 'extreme' channels in group to trim 
D_sel_trial_mask_ref = bml_rereference(cfg,D_sel_trial_mask);


% % % Adding unfiltered channels

%% Adding unfiltered channels
cfg =[];
cfg.channel = setdiff(D.label, D_sel_trial_mask_ref.label);
D_unfilt = ft_selectdata(cfg,D);

cfg = [];
cfg.epoch = epoch;
D_unfilt_trial = bml_redefinetrial(cfg,D_unfilt);

%%% AM got rid of combining unfilt and filt_trial_mask_ref because the time fields do not match
% % % % % % % % % % cfg=[];
% % % % % % % % % % cfg.appenddim = 'chan';
% % % % % % % % % % D_trial_ref = ft_appenddata(cfg,D_sel_filt_trial_mask_ref, D_unfilt_trial);
D_trial_ref = D_sel_trial_mask_ref; 

% % % Saving referenced data

bml_annot_write(epoch,['annot/' SUBJECT '_trial_epoch.txt']);
save([PATH_SAVE_PREPROCESSED filesep SUBJECT '_ft_raw_filt_trial_denoised_ref_criteria_E.mat'],'D_trial_ref','-v7.3');

% % % Quality check - visually inspect the data

cfg=[];
cfg.viewmode = 'vertical';
cfg.blocksize = 8;
cfg.ylim = 'maxmin';
cfg.continuous = 'no';
ft_databrowser(cfg,D_trial_ref);

% % % Quality check - create crosscorrelation matrix
% % % remask with zeros to calculate crosscorrelation

cfg=[];
cfg.remask_nan = true;
cfg.value = 0;
% % % % D_trial_ref_mask0 = bml_mask(cfg,D_trial_ref);
D_trial_ref_mask0 = bml_mask(cfg,D_sel_trial_mask_ref);


% % % do timelock analysis for the raw, high pass filtered (hpf) and rereferenced (ref) versions of the data


cfg=[];
cfg.covariance = 'yes';
cfg.vartrllength = 2;
cfg.trials = 1; %selecting only first trial to assess crosscorrelation
TL_sel=ft_timelockanalysis(cfg,D_sel);

cfg=[];
cfg.covariance = 'yes';
cfg.vartrllength = 2;
cfg.trials = 1; %selecting only first trial to assess crosscorrelation
TL_sel_filt=ft_timelockanalysis(cfg,D_sel);

cfg=[];
cfg.covariance = 'yes';
cfg.vartrllength = 2;
cfg.trials = 1; %selecting only first trial to assess crosscorrelation
TL_ref_mask0=ft_timelockanalysis(cfg,D_trial_ref_mask0);


% % % Plot crosscorrelation matrix for raw, hpf and car objects

f=figure('Position',[0 0 1500 300]);
subplot(1,3,1)
image(corrcov(TL_sel.cov),'CDataMapping','scaled')
caxis([-1 1])
title('raw')
colorbar()

subplot(1,3,2)
image(corrcov(TL_sel_filt.cov),'CDataMapping','scaled')
caxis([-1 1])
title('hpf')
colorbar()

subplot(1,3,3)
image(corrcov(TL_ref_mask0.cov),'CDataMapping','scaled')
caxis([-1 1])
title('ref')
colorbar()

saveas(f,['figures/' SUBJECT '_P09_raw_filt_ref_xcorr_T1.png'])
