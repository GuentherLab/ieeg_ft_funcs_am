%% currently only working with Triplet... need move certain parameters into the calling script to make more generic
%%% do referencing on fieldtrip data that has already been vibration denoised and annotated for artifacts
%%% called by e.g. batch_P08_P09_triplet

function P09_redefine_trial_reref(op)

% % % % % Loading packages
ft_defaults
bml_defaults
format long
set_project_specific_variables() % set paths, load subject-specific variables including electrode and artifact table

%%% these filtering parameters only apply if loading non-denoised data, which will need to be filtered
field_default('op','do_high_pass_filter','yes'); % should a high pass filter be applied
    field_default('op','high_pass_filter_freq',1); %cutoff frequency of high pass filter, if one is to be used
field_default('op','do_bsfilter','yes'); % should a band stop filter be applied
    field_default('op','line_noise_harm_freqs',[60 120 180 240]); % for notch filters for 60hz harmonics, if band stop filter is to be used
field_default('op','out_freq',100); % downsample rate in hz

% filenamename that we will save rereferenced data into
ft_savename = [FT_FILE_PREFIX 'raw_filt_trial_ar-',op.art_crit, '_ref-',op.rereference_method, op.denoise_string, '.mat']; 


%%
% % % Load FieldTrip raw data
fprintf('\n* Doing rereferencing (method = %s; artifact criteria %s) for subject %s...',op.rereference_method, op.art_crit, op.sub)


%% for non denosied data, need to load raw fieldtrip here rather than load trialwise data (which is done for denoised).... feed into redefinetrial

%%%% load raw data...
%%%% for denoised, this must be trialwise 
%%%% for not denoised, use continuous (trialwise raw hasn't necessarily been saved)... filter this data
%%%% .... in either case, the same trial boundaries will be used for cutting fieldtrip data into trial epochs
if op.denoised
    load([FT_FILE_PREFIX 'raw_filt_trial_', op.denoise_string, '.mat'], 'D');
    D_filt = D; clear D    % denoised data is already filtered
elseif ~op.denoised 
    load([FT_FILE_PREFIX 'raw.mat'], 'D', 'loaded_epoch');

   % Applying High Pass Filter and remove 60hz line noise
    cfg=[];
    cfg.hpfilter=op.do_high_pass_filter;
        cfg.hpfreq=op.high_pass_filter_freq;
    cfg.hpfilttype='but';
    cfg.hpfiltord=5;
    cfg.hpfiltdir='twopass';




%     cfg.bsfilter = op.do_bsfilter;
%         cfg.bsfreq= [op.line_noise_harm_freqs-1; op.line_noise_harm_freqs+1]'; % notch filters for 60hz harmonics
%     

    cfg.bsfilter = 'no'; 
    cfg.dftfilter='yes';
    cfg.dftfreq           = [60 120 180 240 300 360 420 480];
    cfg.dftbandwidth      = [ 1   1   1   1   1   1   1   1];
    cfg.dftneighbourwidth = [ 2   2   2   2   2   2   2   2];


    cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'};
    D_filt = ft_preprocessing(cfg,D);
end

% % % Selecting electrodes and remasking with zeros instead of NaNs
% % % Some of FieldTrip's functions don't work with NaNs, so we are going to temporarily replace NaNs with zeros to avoid issues. 

cfg1=[];
cfg1.channel={'ecog_*','macro_*','micro_*','dbs_*'};
D_sel_filt = ft_selectdata(cfg1,D_filt);

%% Redefining trials
cfg1 = [];
cfg1.epoch = epoch; % epoch should be set in set_project_specific_variables - determines trial boundaries
D_sel_filt_trial = bml_redefinetrial(cfg1,D_sel_filt); % D was already in trial format, but we are redefining trial timing


% % % Masking artifacts and empty electrodes with NaNs for re-referencing

%combining artifacts with empty_electrode table
empty_electrode_file = [PATH_ANNOT, filesep, op.sub '_empty_electrode.txt']; 
if exist('empty_electrode_file', 'file')
    empty_electrode = bml_annot_read();
    cfg1=[];
    cfg1.groupby='label';
    artifact_and_empty = bml_annot_union(cfg1, artifact, empty_electrode);
else
    artifact_and_empty = artifact;
end

    %masking and empty_electrodes with NaNs
cfg1=[];
cfg1.annot=artifact_and_empty;
cfg1.label_colname = 'label';
cfg1.complete_trial = true; %masks entire trials
cfg1.value=NaN;
D_sel_filt_trial_mask = bml_mask(cfg1,D_sel_filt_trial);

ieeg_rows = cellfun(@(x)any(strcmp(x,{'ecog','ECOG','dbs','DBS','macro','SEEG'})), electrodes.type);
el_ieeg = electrodes(ieeg_rows,:);

if strcmp(op.rereference_method,'none') % no referencing
    D_sel_filt_trial_mask_ref = D_sel_filt_trial_mask;
elseif ~strcmp(op.rereference_method,'none')
    cfg=[];
    cfg.label = el_ieeg.name;
    cfg.group = el_ieeg.connector;
    cfg.method = op.rereference_method; % 
    cfg.percent = 50; %percentage of 'extreme' channels in group to trim 
    D_sel_filt_trial_mask_ref = bml_rereference(cfg,D_sel_filt_trial_mask);
end


% % % Adding unfiltered channels

%% Adding unfiltered channels
%%% AM got rid of combining unfilt and filt_trial_mask_ref because the time fields do not match
% % % % % % % % % % cfg1 =[];
% % % % % % % % % % cfg1.channel = setdiff(D.label, D_sel_trial_mask_ref.label);
% % % % % % % % % % D_unfilt = ft_selectdata(cfg1,D);
% % % % % % % % % % 
% % % % % % % % % % cfg1 = [];
% % % % % % % % % % cfg1.epoch = epoch;
% % % % % % % % % % D_unfilt_trial = bml_redefinetrial(cfg1,D_unfilt);


% % % % % % % % % % cfg=[];
% % % % % % % % % % cfg.appenddim = 'chan';
% % % % % % % % % % D_trial_ref = ft_appenddata(cfg,D_sel_filt_trial_mask_ref, D_unfilt_trial);
D_trial_ref = D_sel_filt_trial_mask_ref; 

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
D_trial_ref_mask0 = bml_mask(cfg1,D_trial_ref);


% % % do timelock analysis for the raw, high pass filtered (hpf) and rereferenced (ref) versions of the data


cfg1=[];
cfg1.covariance = 'yes';
cfg1.vartrllength = 2;
cfg1.trials = 1; %selecting only first trial to assess crosscorrelation
TL_sel=ft_timelockanalysis(cfg1,D_sel_filt);

cfg1=[];
cfg1.covariance = 'yes';
cfg1.vartrllength = 2;
cfg1.trials = 1; %selecting only first trial to assess crosscorrelation
TL_sel_filt=ft_timelockanalysis(cfg1,D_sel_filt);

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