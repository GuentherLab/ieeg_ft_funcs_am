%%% do offline rereferencing of raw [HPF and line-noise filtered] signal after applying artifact masks
%%% ... the function P08A09_detect_artifact_not_denoised must have been run first
% protocol P09 for Pitt data; not yet assigned a protocol number for MGH data

function P09_redefine_trial_common_avg_ref_not_denoised(op)

% Loading packages
ft_defaults
bml_defaults
format long

%% Defining paths, loading artifact parameters
vardefault('op',struct); % initialize options if not present
field_default('op','sub','DM1005')
field_default('op','art_crit','E'); % 'E' = 70-250hz high gamma; 'F' = beta; 'G' = other Rohan criterion? 
field_default('op','do_high_pass_filter','yes'); % should a high pass filter be applied
    field_default('op','high_pass_filter_freq',1); %cutoff frequency of high pass filter, if one is to be used
field_default('op','rereference_method','CTAR'); 
field_default('op','out_freq',100); % downsample rate in hz for high gamma traces

% use the following trial duration if one is not specified in the trials annot table
%%% ..... probably not specified due to a missing keypress at the end of a session
% if not the last trial, use default_iti for trial spacing
field_default('op','default_trialdur_max_if_empty',15);  % sec
field_default('op','default_iti_if_empty',0.5); % sec

set_project_specific_variables(); % set paths etc. based on data collection site, load timing and electrode data

% filenamename that we will save rereferenced data into
ft_savename = [FT_FILE_PREFIX 'raw-filt-trial_ar-',op.art_crit, '_ref-',op.rereference_method,  '_not-denoised.mat']; 

% handle missing trial durations
for itrial = 1:height(epoch)
    if isnan(epoch.duration(itrial))
        if ~[itrial == height(epoch)]
            epoch.ends(itrial)= min([epoch.starts(itrial) + op.default_trialdur_max_if_empty, epoch.ends(itrial+1) - op.default_iti_if_empty]);
        elseif itrial == height(epoch) % last trial
            epoch.ends(itrial)= epoch.starts(itrial) + op.default_trialdur_max_if_empty;
        end
        epoch.duration(itrial) = epoch.ends(itrial) - epoch.starts(itrial); 
    end
end

% % % Load FieldTrip raw data
load(FT_RAW_FILENAME,'D','loaded_epoch');

% % % Adjusting length of sessions for notch filter to work
D_annot = bml_raw2annot(D);
Fl=[60 120 180 240 300 360 420 480];
D_annot.nSamples2 = round(floor(D_annot.nSamples .* Fl(1)./D_annot.Fs) .* D_annot.Fs./Fl(1));
D_annot.nSamples2 = round(floor(D_annot.nSamples2 .* Fl(1)./D_annot.Fs) .* D_annot.Fs./Fl(1));
D_annot.nSamples2 = D_annot.nSamples2(:,1);
D_annot.ends = D_annot.starts + D_annot.nSamples2 ./ D_annot.Fs;

cfg=[];
cfg.epoch=D_annot;
D1= bml_redefinetrial(cfg,D);


% % % Selecting electrodes and remasking with zeros instead of NaNs

% % % Some of FieldTrip's functions don't work with NaNs, so we are going to temporarily replace NaNs with zeros to avoid issues. 

cfg=[];
cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'};
% % %         cfg.trials = logical([1 1 1 1 0]);
D_sel = ft_selectdata(cfg,D1);

cfg=[];
cfg.remask_nan = true;
cfg.value = 0;
D_sel = bml_mask(cfg, D_sel);

% % % Applying high pass filter and line noise removal filter
% AM note: this step uses dft filter instead of bandstop filter (like artifact rejection uses).... ask Alan why the difference
cfg=[];
cfg.hpfilter=op.do_high_pass_filter;
    cfg.hpfreq=op.high_pass_filter_freq;
cfg.hpfilttype='but';
cfg.hpfiltord=5;
cfg.hpfiltdir='twopass';
cfg.dftfilter='yes';
%%%% this interpolation option is causing an error with 3012, artifact criterion E
% cfg.dftreplace='neighbour';  %using spectrum interpolation method Mewett et al 2004
cfg.dftfreq           = [60 120 180 240 300 360 420 480];
cfg.dftbandwidth      = [ 1   1   1   1   1   1   1   1];
cfg.dftneighbourwidth = [ 2   2   2   2   2   2   2   2];
D_sel_filt = ft_preprocessing(cfg,D_sel);

% % % Redefining trials
% for dbs-seq/smsl, we will use experimenter keypress for trial start/end times
%%%%% this means no trial overlap, but generally a large time buffer before cue onset and after speech offset
cfg = [];
cfg.epoch = epoch;
D_sel_filt_trial = bml_redefinetrial(cfg,D_sel_filt);


% % % Masking artifacts with NaNs for re-referencing

%masking artifacts and empty_electrodes with NaNs
cfg=[];
cfg.annot=artifact;
cfg.label_colname = 'label';
cfg.complete_trial = true; %masks entire trials
cfg.value=NaN;
D_sel_filt_trial_mask = bml_mask(cfg,D_sel_filt_trial);

% % % Common trimmed average reference per connector groups

el_ecog = electrodes(electrodes.type=="ECOG" | electrodes.type=="ecog",:);

if strcmp(op.rereference_method,'none') % no referencing
    D_sel_filt_trial_mask_ref = D_sel_filt_trial_mask;
elseif ~strcmp(op.rereference_method,'none')
    cfg=[];
    cfg.label = el_ecog.name;
    cfg.group = el_ecog.connector;
    cfg.method = op.rereference_method; % 
    cfg.percent = 50; %percentage of 'extreme' channels in group to trim 
    D_sel_filt_trial_mask_ref = bml_rereference(cfg,D_sel_filt_trial_mask);
end


% % % Adding unfiltered channels

%% Adding unfiltered channels
cfg =[];
cfg.channel = setdiff(D.label, D_sel_filt_trial_mask_ref.label);
D_unfilt = ft_selectdata(cfg,D);

cfg = [];
cfg.epoch = epoch;
D_unfilt_trial = bml_redefinetrial(cfg,D_unfilt);

%%% for Triplet task, AM got rid of combining unfilt and filt_trial_mask_ref because the time fields do not match
% % % % % % % % % % cfg=[];
% % % % % % % % % % cfg.appenddim = 'chan';
% % % % % % % % % % D_trial_ref = ft_appenddata(cfg,D_sel_filt_trial_mask_ref, D_unfilt_trial);
D_trial_ref = D_sel_filt_trial_mask_ref; 

% % % Saving referenced data
% this section used to saved an annot table file called ['annot/' SUBJECT '_trial_epoch.txt'] from epoch variable... AM removed it 2024/02/05 because it appeared unnecessary and confusing
save(ft_savename,'D_trial_ref','-v7.3');

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
D_trial_ref_mask0 = bml_mask(cfg,D_sel_filt_trial_mask_ref);


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
TL_sel_filt=ft_timelockanalysis(cfg,D_sel_filt);

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

saveas(f,[PATH_FIGURES filesep 'sub-' SUBJECT '_P09_raw_filt_ref_xcorr_T1_not_denoised.png'])
