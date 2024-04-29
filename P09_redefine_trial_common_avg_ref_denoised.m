%% currently only working with Triplet... need move certain parameters into the calling script to make more generic
%%% do referencing on fieldtrip data that has already been vibraiton denoised and annotated for artifacts
%%% called by e.g. batch_P08_P09_triplet

function P09_redefine_trial_common_avg_ref_denoised(op)

% % % % % Loading packages
ft_defaults
bml_defaults
format long
set_project_specific_variables() % set paths, load subject-specific variables including electrode and artifact table

% filenamename that we will save rereferenced data into
ft_savename = [FT_FILE_PREFIX 'raw_filt_trial_ar-',op.art_crit, '_ref-',op.rereference_method, op.denoise_string, '.mat']; 

% padding time to include before/after trial onset/offset
field_default('op','prebuffer_sec',1.5); 
field_default('op','postbuffer_sec',2); 

% % % Load annotation tables
empty_electrode = bml_annot_read([PATH_ANNOT, filesep, op.sub '_empty_electrode.txt']);
cue_presentation = bml_annot_read([PATH_ANNOT, filesep, op.sub '_cue_presentation.txt']);

%%
% % % Load FieldTrip raw data
fprintf('\n* Doing rereferencing (method=%s; artifact criteria %s) for subject %s...',op.rereference_method, op.art_crit, op.sub)
load([FT_FILE_PREFIX 'raw_filt_trial_denoised.mat']);

% % % Selecting electrodes and remasking with zeros instead of NaNs
% % % Some of FieldTrip's functions don't work with NaNs, so we are going to temporary replace NaNs with zeros to avoid issues. 

cfg1=[];
cfg1.channel={'ecog_*','macro_*','micro_*','dbs_*'};
% % %         cfg.trials = logical([1 1 1 1 0]);
D_sel = ft_selectdata(cfg1,D);


%% Redefining trials

% % % If the cue_presentation table has other column names, modify the following lines accordingly. 
% % % Use generous buffers around the cue and speech production, even if this means that consecutive trials overlap. 

epoch = cue_presentation(:,{'stim1_onset','ends','session_id','trial_id'});
epoch.starts = epoch.stim1_onset - op.prebuffer_sec;
epoch.ends = epoch.ends + op.postbuffer_sec;
epoch = bml_annot_table(epoch);

cfg1 = [];
cfg1.epoch = epoch;
D_sel_trial = bml_redefinetrial(cfg1,D_sel); % D was already in trial format, but we are redefining trial timing


% % % Masking artifacts and empty electrodes with NaNs for re-referencing

%combining artifacts with empty_electrode table
cfg1=[];
cfg1.groupby='label';
artifact_empty = bml_annot_union(cfg1, artifact, empty_electrode);

%masking artifacts and empty_electrodes with NaNs
cfg1=[];
cfg1.annot=artifact_empty;
cfg1.label_colname = 'label';
cfg1.complete_trial = true; %masks entire trials
cfg1.value=NaN;
D_sel_trial_mask = bml_mask(cfg1,D_sel_trial);

% % % Common trimmed average reference per connector groups

el_ecog = electrodes(electrodes.type=="ecog",:);

cfg1=[];
cfg1.label = el_ecog.electrode;
cfg1.group = el_ecog.connector;
cfg1.method = 'CTAR'; %using trimmed average referencing
cfg1.percent = 50; %percentage of 'extreme' channels in group to trim 
D_sel_trial_mask_ref = bml_rereference(cfg1,D_sel_trial_mask);


% % % Adding unfiltered channels

%% Adding unfiltered channels
cfg1 =[];
cfg1.channel = setdiff(D.label, D_sel_trial_mask_ref.label);
D_unfilt = ft_selectdata(cfg1,D);

cfg1 = [];
cfg1.epoch = epoch;
D_unfilt_trial = bml_redefinetrial(cfg1,D_unfilt);

%%% AM got rid of combining unfilt and filt_trial_mask_ref because the time fields do not match
% % % % % % % % % % cfg=[];
% % % % % % % % % % cfg.appenddim = 'chan';
% % % % % % % % % % D_trial_ref = ft_appenddata(cfg,D_sel_filt_trial_mask_ref, D_unfilt_trial);
D_trial_ref = D_sel_trial_mask_ref; 

% % % Saving referenced data

bml_annot_write(epoch,[PATH_ANNOT, filesep, op.sub '_trial_epoch.txt']);
save(ft_savename,'D_trial_ref','-v7.3');

% % % Quality check - visually inspect the data

cfg1=[];
cfg1.viewmode = 'vertical';
cfg1.blocksize = 8;
cfg1.ylim = 'maxmin';
cfg1.continuous = 'no';
ft_databrowser(cfg1,D_trial_ref);

% % % Quality check - create crosscorrelation matrix
% % % remask with zeros to calculate crosscorrelation

cfg1=[];
cfg1.remask_nan = true;
cfg1.value = 0;
% % % % D_trial_ref_mask0 = bml_mask(cfg,D_trial_ref);
D_trial_ref_mask0 = bml_mask(cfg1,D_sel_trial_mask_ref);


% % % do timelock analysis for the raw, high pass filtered (hpf) and rereferenced (ref) versions of the data


cfg1=[];
cfg1.covariance = 'yes';
cfg1.vartrllength = 2;
cfg1.trials = 1; %selecting only first trial to assess crosscorrelation
TL_sel=ft_timelockanalysis(cfg1,D_sel);

cfg1=[];
cfg1.covariance = 'yes';
cfg1.vartrllength = 2;
cfg1.trials = 1; %selecting only first trial to assess crosscorrelation
TL_sel_filt=ft_timelockanalysis(cfg1,D_sel);

cfg1=[];
cfg1.covariance = 'yes';
cfg1.vartrllength = 2;
cfg1.trials = 1; %selecting only first trial to assess crosscorrelation
TL_ref_mask0=ft_timelockanalysis(cfg1,D_trial_ref_mask0);


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

saveas(f,[PATH_FIGURES, filesep, op.sub '_P09_raw_filt_ref_xcorr_T1.png'])
