%% currently only working with Triplet


% AM 

% % % Do a time-frequency analysis
% % % 
% % % Save the working script as <subject_id>/Preprocessed Data/Sync/<subject_id>_P10_timefreq_analysis_<date_yyyymmdd>.m. 
% % % 
% % % Make sure FieldTrip and BML folders are in the path
% % % Load defaults (adds additional folders to Matlab’s path) and define paths to files

ft_defaults
bml_defaults
format long

SUBJECT='DBS3012';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
ARTIFACT_CRIT = 'E';
cd(PATH_SYNC)

% % % Load annotation tables

coding = bml_annot_read(['annot/' SUBJECT '_coding.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
artifact = bml_annot_read(['annot/' SUBJECT '_artifact_criteria_', ARTIFACT_CRIT, '.txt']);
empty_electrode = bml_annot_read(['annot/' SUBJECT '_empty_electrode.txt']);

% % % Load FieldTrip common average referenced data

load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep...
    SUBJECT '_ft_raw_filt_trial_ar_ref_criteria_', ARTIFACT_CRIT, '.mat']);

% % % Defining table to epoch the data

epoch = table(); %creating empty table
epoch.syl1_onset = coding.syl1_onset; %copying onset of syllable 1 from coding table
epoch.trial_id = coding.trial_id; 
epoch.session_id = coding.session_id; 
epoch.starts = coding.syl1_onset - 2; %creating starts
epoch.ends = coding.syl1_onset + 3; %creating ends
epoch = bml_annot_table(epoch); %creating id and duration columns
epoch = epoch(~ismissing(epoch.starts),:); %removing rows with NaNs


% % % Checking the epoching by the audio envelopes

cfg=[];
cfg.channel = 'envaudio_*';
D_trial_ref_audio = ft_selectdata(cfg,D_trial_ref);

cfg=[];
cfg.epoch = epoch;
cfg.timelock = 'syl1_onset';
cfg.timesnap = true;
D_trial_ref_audio = bml_redefinetrial(cfg,D_trial_ref_audio);

cfg=[];
cfg.ylim = 'maxmin';
% ft_databrowser(cfg,D_trial_ref_audio);

% % % Select channels on which to perform time frequency analysis
% % % 	In this example we are selecting the first ecog strip only

cfg=[];
cfg.channel = {'ecog_*','dbs_*','micro_*','macro_*','audio_*'};
D_trial_ref_sel = ft_selectdata(cfg,D_trial_ref);

%% Redefining trials based on new epochs. 
% Timelocking to first syllable onset. 
% "Timesnap" the time vector to avoid numerical issues. 
% For this dataset, this means rounding the time vector to the closest millisecond.  

% redefining trials and timelocking
cfg=[];
cfg.epoch = epoch;
cfg.timelock = 'syl1_onset';
cfg.timesnap = true;
D_trial_ref_sel_epoch = bml_redefinetrial(cfg,D_trial_ref_sel);

%% Mask entire trials for electrodes, when artifacts were detected.

% remasking NaNs with zeros
cfg=[];
cfg.value=0;
cfg.remask_nan = true;
cfg.complete_trial = true;
D_trial_ref_sel_epoch = bml_mask(cfg,D_trial_ref_sel_epoch);

% % % Redefine trials using the epoch table and the masked data object

% % % % % % % cfg=[];
% % % % % % % cfg.epoch = epoch;
% % % % % % % cfg.t0 = 'onset_syl1';
% % % % % % % D_ref_ecog1_mask0_epoch = bml_redefinetrial(cfg,D_ref_ecog1_mask0);

%%  Perform the time frequency analysis
% make take tens of minutes to run

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'wavelet';
cfg.foi          = 2:2:250;   % analysis 2 to 250 Hz in steps of 2 Hz 
cfg.toi          = -1.5:0.01:2.5; %from -1.5 to 2.5 sec in steps of 10 ms
cfg.feedback     = 'no';
cfg.pad          = 'nextpow2';
TF_trial_ref_sel_epoch = ft_freqanalysis(cfg, D_trial_ref_sel_epoch); 

%% Load the ECoG strip layout and plot
cfg=[];
cfg.layout = [PATH_SYNC filesep 'annot' filesep SUBJECT '_ecog2.lay'];
cfg.skipscale = 'yes';
cfg.outline = 'no';
cfg.mask = 'no';
lay = ft_prepare_layout(cfg);
lay.width = lay.width*2.95;
lay.height = lay.height*2.7;

cfg = [];
cfg.baseline     = [-2 -1.5]; 
cfg.baselinetype = 'db'; 
cfg.xlim         = [-0.5 2.5];
cfg.ylim         = [0 250];
cfg.zlim         = [-1 6];	        
cfg.showlabels   = 'no';	
cfg.layout       = lay;
cfg.colorbar     = 'yes';
cfg.showscale    = 'no';
cfg.showcomment  = 'no';
f=figure('Position',[0 0 1000 2000]);
ft_multiplotTFR(cfg, TF_trial_ref_sel_epoch);
saveas(f,['figures/' SUBJECT '_P10_timefreq_analysis_db_ecog1_SMC.png'])


%% Load the DBS lead layout and plot

cfg=[];
cfg.layout = [PATH_SYNC '\annot\' SUBJECT '_dbs_Medtronic_4.lay'];
cfg.skipscale = 'yes';
cfg.outline = 'no';
cfg.mask = 'no';
lay = ft_prepare_layout(cfg);
lay.width = lay.width*3;
lay.height = lay.height*0.9;

cfg = [];
cfg.baseline     = [-2 -1.5]; 
cfg.baselinetype = 'db'; 
cfg.xlim         = [-0.5 2.5];
cfg.ylim         = [0 250];
cfg.zlim         = [-1 4];	        
cfg.showlabels   = 'no';	
cfg.layout       = lay;
cfg.colorbar     = 'yes';
cfg.showscale    = 'no';
cfg.channel      = 'dbs_*';
cfg.showcomment  = 'no';
f=figure('Position',[0 0 1*300 4*100]);
ft_multiplotTFR(cfg, TF_trial_ref_sel_epoch);
saveas(f,['figures/' SUBJECT '_P10_timefreq_analysis_db_dbs_STN.png'])

%% Load macro and micro layouts and plot
% Macro/Micro electrodes
%plotting aggregated sessions
cfg=[];
cfg.layout = [PATH_SYNC '\annot\macro_micro.lay'];
cfg.skipscale = 'yes';
cfg.outline = 'no';
cfg.mask = 'no';
lay = ft_prepare_layout(cfg);
lay.pos(end-1,1)=-5;
lay.width = lay.width*2.95;
lay.height = lay.height*2.7;

cfg = [];
cfg.baseline     = [-2 1.5]; 
cfg.baselinetype = 'db'; 
cfg.xlim         = [-0.5 2.5];
cfg.ylim         = [0 250];      
cfg.showlabels   = 'no';
cfg.layout       = lay;
cfg.colorbar     = 'yes';
cfg.showscale    = 'no';
cfg.showcomment  = 'yes';
cfg.channel      = {'macro_*','micro_*'};

