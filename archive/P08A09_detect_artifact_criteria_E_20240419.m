%% currently only working with Triplet

function P08A09_detect_artifact_criteria_E(SUBJECT, param)

%%% Detect artifacts.... for use with vibration-denoised data
% protocol P08 in Pitt data, protocol A09 in MGH data
%%%% before running this function, run P08A09_highgamma_from_denoised to extract trialwise high gamma
%%%%.... these trials will be concatenated in this script before performing artifact detection

% creates an artifact annotation table

% % % % CRITERIA E paramaters valus
CRITERIA = 'E'; %identifier for the criteria implemented in this script

% % % % % load packages
ft_defaults
bml_defaults
format long

% % % % % % Defining paths
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_PROTOCOL = 'Z:\DBS\Batch\P08_artifact_criteria_E';
PATH_ANNOT = [PATH_SYNC '/annot']; 

cd(PATH_SYNC)

session= bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);

% DBS3031 is missing DBS3031_trial_epoch.txt in annot folder; instead only has DBS3031_trial_epoch_criteria_D.txt
if ~strcmp(SUBJECT, 'DBS3031') 
    trial_epoch = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch.txt']);
elseif strcmp(SUBJECT, 'DBS3031') 
    trial_epoch = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch_criteria_D.txt']);
end

%% Loading FieldTrip data 
load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_hg_trial_denoised.mat'],'D_hg_trial');
% load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_hg_trial_denoised.mat'],'D_wavpow_trial');


%   cut out of trialtable all rows for which the vibration-denoised trial is missing
cfg = []; 
cfg.plot_times = 0; 
cfg.trials = trial_epoch; 
[trials, trials_ft] = P08_correct_fieldtrip_trialtable_discrepancies(cfg,D_hg_trial); 

%% concatenate non-overlaping parts of trials, sorting by run
session_ids = unique(trials_ft.session_id);
nses = length(session_ids); 
ntrials = length(D_hg_trial.trial);
nchans = length(D_hg_trial.label); 
tt_overlapping = D_hg_trial.time;
sampint = 1/D_hg_trial.fsample;

    %%%%%% get new trial times that will not overlap with each other
trialtimes_no_overlap = table(nan(ntrials,1),nan(ntrials,1),'VariableNames',{'starts','ends'}); 
interp_times = [];
for itrial = 1:ntrials-1
    next_trial_start = tt_overlapping{itrial+1}(1);
    nonoverlap_timepoints = tt_overlapping{itrial}(tt_overlapping{itrial} < next_trial_start); 

    %%% if this trial has any timepoints non-overlapping with the subsequent trials, mark that nonoveralpping window
    %%% .... if it's entirely contained with the subsequent trial, leave its window as nans
    if ~isempty(nonoverlap_timepoints)
  
        % trial onsets stay the same
        trialtimes_no_overlap.starts(itrial) = tt_overlapping{itrial}(1); 
    
        % for trial offsets, take the latest timepoint of this trial which is....
        % ... before the first timepoint of the subsequent trial
        trialtimes_no_overlap.ends(itrial) = max(nonoverlap_timepoints); 
        
        % look for time discontinuities between trials within a session
        if trials_ft.session_id(itrial) == trials_ft.session_id(itrial+1) && ... % if trials are in the same session
                next_trial_start - trialtimes_no_overlap.ends(itrial) > 2*sampint % ...and intertrial gap exceeds 2 timesteps
            interp_times_this_trial = trialtimes_no_overlap.ends(itrial) + sampint : sampint : next_trial_start;
            interp_times = [interp_times, interp_times_this_trial]; 
        end
    end
end
trialtimes_no_overlap.starts(ntrials) = tt_overlapping{end}(1); % last trial times are unchanged
trialtimes_no_overlap.ends(ntrials) = tt_overlapping{end}(end); % last trial times are unchanged
trialtimes_no_overlap.duration = trialtimes_no_overlap.ends - trialtimes_no_overlap.starts;

% re-cut HG trialwise data with non-overlapping trial times
%%% there should not be gaps more than 1 timestep between redefinied trials....
%%%     ... (except possibly if there was a long pause betwteen trials within a session)
cfg = [];
cfg.epoch = trialtimes_no_overlap;
D_hg_trial_no_overlap = bml_redefinetrial(cfg,D_hg_trial);


%%%%% concatenate trials within each session to recreate continuous session-wise data
D_hg = struct;
D_hg.label = D_hg_trial.label;
D_hg.trial = cell(1,nses);
D_hg.time = cell(1,nses); 
D_hg.fsample = D_hg_trial.fsample;

%%% get session start/end times from trial times
cfg = [];
cfg.plot_times = 0; 
cfg.trials = trials; 
[~, trials_ft_no_overlap] = P08_correct_fieldtrip_trialtable_discrepancies(cfg,D_hg_trial_no_overlap); % remake trials_ft in case some trials got removed
sestab = table(session_ids, nan(nses,1), nan(nses,1), nan(nses,1), 'VariableNames', {'session_id', 'starts', 'ends', 'duration'}); 
for ises = 1:nses
    trialmatch = find(trials_ft_no_overlap.session_id == session_ids(ises)); % row indices
    sestab.starts(ises) = D_hg_trial_no_overlap.time{trialmatch(1)}(1); % start of the first trial of this session
    sestab.ends(ises) = D_hg_trial_no_overlap.time{trialmatch(end)}(end); % end of the last trial of this session
    times_pre_interp  = cell2mat(D_hg_trial_no_overlap.time(trialmatch));
    % times to interpolate this session
    times_interp_session = interp_times( interp_times > sestab.starts(ises) & interp_times < sestab.ends(ises));
    times_with_interp = [times_pre_interp, times_interp_session]; 
    [D_hg.time{ises}, timeorder] = sort(times_with_interp); % sort timepoints
    ephys_pre_interp = cell2mat(D_hg_trial_no_overlap.trial(trialmatch)); 
    ephys_interp_session = nan(nchans, length(times_interp_session)); 
    ephys_with_interp = [ephys_pre_interp, ephys_interp_session]; % add interpolated points as NaNs
    D_hg.trial{ises} = ephys_with_interp(:,timeorder); % reorder ephys data according to time
end
sestab.duration = sestab.ends - sestab.starts; 



%% loading electrode type band table

% working in protocol folder
cd(PATH_PROTOCOL)

if ~exist('el_band','var')
  param = readtable('artifact_E_params.txt');
  param_default = param(param.subject == "default",:);
  param_subject = param(strcmp(param.subject,SUBJECT),:);
  if ~isempty(param_subject)
    param = bml_annot_rowbind(param_default(~ismember(param_default.name,param_subject.name),:),param_subject);
  end
end

%%%% we are working with high-gamma power data (not raw), so don't do highpass or line noise filters

%% Artifact rejection
% iterating over bands and electrode types
artifact = table();
f = figure();
for idx = 1:height(param)
  
  fprintf('doing %s %s \n',SUBJECT,param.name{idx});
  
  el_type = strip(param.electrode_type{idx});
  env_mult_factor =  param.env_mult_factor(idx);
  pname = strip(param.name{idx});
  
  ENVELOPE_BIN_SIZE_SECONDS = param.env_bin_size(idx); %envelope bin size in seconds
  THRESHOLD_STD_FACTORS = [param.th_factor_std_low(idx), param.th_factor_std_high(idx)]; %factors to determine detection thresholds 
  THRESHOLD_FIX = [param.th_fix_min(idx), param.th_fix_max(idx)]; %fix thresholds to filter data before applying robust estimates
  CONSOLIDATION_TIME_TOLERANCE = param.th_consolidation(idx); %min time allowed between consecutive artifacts
  ELECTRODE_COVERAGE_THRESHOLD = param.th_frac_coverage(idx); %max allowed fraction of time with artifacts
  CONNECTOR_THRESHOLD = [param.th_conn_low(idx), param.th_conn_high(idx)]; %detection threshold for number of electrodes in a connector  
  
  %selecting ECoG channels for artifact rejection
  cfg=[];
  cfg.channel = [el_type,'_*'];
  D_hg_eltype = ft_selectdata(cfg,D_hg);

  if isempty(D_hg_eltype.label)
    %channel type not available
    continue
  end
  
%   % visually inspect signal
%   cfg=[];
%   cfg.viewmode = 'vertical';
%   cfg.blocksize = 30;
%   cfg.ylim = 'maxmin';4
%   cfg.continuous = 'yes';
%   cfg.channel = {'ecog_11*', 'audio_*'};
%   ft_databrowser(cfg,D);

%   cfg=[];
%   cfg.viewmode = 'vertical';
%   cfg.blocksize = 30;
%   cfg.ylim = 'maxmin';
%   cfg.continuous = 'yes';
%   ft_databrowser(cfg,D_hg_ecog_env);

  cfg=[];
  cfg.freq = 1/ENVELOPE_BIN_SIZE_SECONDS;
  D_hg_eltype_env = bml_envelope_binabs(cfg,D_hg_eltype);

  %calculating log10 transform (envelopes have log normal distributions)
  D_hg_eltype_env_log10 = bml_apply(@(x) env_mult_factor .* log10(x),D_hg_eltype_env);
  
  cfg=[];
  cfg.remask_inf=true;
  cfg.value=NaN;
  D_hg_eltype_env_log10 = bml_mask(cfg,D_hg_eltype_env_log10);
  
  %calculating distribution robust statistics. 
  THRESHOLD = nan(nses,2);
  max_v=nan(1,nses);
  min_v=nan(1,nses);
  for i=1:nses
    v = reshape(D_hg_eltype_env_log10.trial{i},1,[]);
    v1 = v((v>THRESHOLD_FIX(1)) & (v<THRESHOLD_FIX(2)));
    m = median(v1);
    std = bml_robust_std(v1);
    if ~isempty(v1)
      max_v(i)=max(v);
      min_v(i)=min(v);
      THRESHOLD(i,:) = m + THRESHOLD_STD_FACTORS.*std;
    end
  end

  %plotting histogram to asses threshold levels
  clf(f); set(f,'Position',[0 0 600 600]);
  for i=1:nses
    subplot(ceil(nses/2),2,i)
    hold on;
    h=histogram(D_hg_eltype_env_log10.trial{i},linspace(min(min_v),max(max_v),61),...
      'FaceAlpha',0.1,'EdgeAlpha',1);
    maxBinCount = max(h.BinCounts);
    plot([THRESHOLD(i,1),THRESHOLD(i,1)],[0,maxBinCount .* 1.1]);
    plot([THRESHOLD(i,2),THRESHOLD(i,2)],[0,maxBinCount .* 1.1]);
    %set(gca,'YScale','log')
    title(['session ' num2str(i)]);
  end
  saveas(f,['figures/' SUBJECT '_' pname '_artifact_env_log10_hist.png'])

  %detecting segments of time for each channel above threshold
  artifact_eltype_1 = table();
  for ises=1:nses
    cfg=[];
    cfg.threshold = THRESHOLD(ises,:);
    cfg.trials = ises;
    artifact_eltype_1 = bml_annot_rowbind(artifact_eltype_1, bml_annot_detect(cfg,D_hg_eltype_env_log10));
  end

  if isempty(artifact_eltype_1)
    continue
  end
  
  %making figure with random snippets of detected artifacts
  cfg=[];
  cfg.n = 1;
  cfg.groupby  = 'label';
  artifact_eltype_1_sample = bml_annot_sample(cfg, artifact_eltype_1);
  artifact_eltype_1_sample = bml_annot_extend(artifact_eltype_1_sample,2,2);

  cfg=[];
  cfg.n = 60;
  artifact_eltype_1_sample = bml_annot_sample(cfg, artifact_eltype_1_sample);
  
  cfg=[];
  cfg.epoch = artifact_eltype_1_sample;
  [D_hg_eltype_sample, epoch_hg_eltype_sample] = bml_redefinetrial(cfg,D_hg_eltype);
 
  D_p = D_hg_eltype_sample;
  E_p = epoch_hg_eltype_sample;
  nx=10; ny=floor(numel(D_p.trial)/nx);
  if ny==0
    ny=1; nx=numel(D_p.trial);
  end
  clf(f); set(f,'Position',[0 0 nx*200 ny*200]);
  for i=1:ny
      for j=1:nx
          pidx = (i-1)*nx+j;
          l = E_p.label(pidx);
          l_idx = bml_getidx(l,D_p.label);
          subplot(ny,nx,pidx);
          plot(D_p.time{pidx},D_p.trial{pidx}(l_idx,:));
          title(E_p.label(pidx));
      end
  end
  saveas(f,['figures/' SUBJECT '_' pname '_artifact_snippets.png'])

  
  %consolidating annotations with CONSOLIDATION_TIME_TOLERANCE margin of overlap
  cfg=[];
  cfg.criterion = @(x) (x.starts(end) - max(x.ends(1:(end-1))) < CONSOLIDATION_TIME_TOLERANCE);
  cfg.groupby = 'label';
  artifact_eltype_2 = bml_annot_consolidate(cfg,artifact_eltype_1);

%   %creating ft_raw from annotations for visualization
%   cfg=[];
%   cfg.template = D_hg_eltype_env_log10;
%   cfg.annot_label_colname='label';
%   artifact_eltype_3_raw = bml_annot2raw(cfg,artifact_eltype_2);
% 
%   %raster plot of artifacts for session 1
%   f=figure();
%   bml_plot_raster(artifact_eltype_3_raw)

  %check if excluded segments are correct
%   cfg=[];
%   cfg.label_colname = 'label';
%   cfg.annot = artifact_eltype_2;
%   cfg.value = NaN;
%   D_hg_eltype_mask = bml_mask(cfg, D_hg_eltype);
% 
%   cfg=[];
%   cfg.viewmode = 'vertical';
%   cfg.blocksize = 30;
%   cfg.ylim = 'maxmin';
%   cfg.continuous = 'yes';
%   ft_databrowser(cfg,D_hg_eltype_mask);

  %% rejecting faulty channels 

  %decide which artifacts to include. Usually just ECoG artifacts
  %artifact = bml_annot_rowbind(artifact_ecog_3,artifact_macro_3,artifact_dbs_3);
  artifact_1 = artifact_eltype_2;


  cfg = [];
  cfg.groupby = 'label';
  artifact_1_session_cvg = bml_annot_coverage(cfg,artifact_1,session);

  %histogram(artifact_1_session_cvg.coverage,linspace(0,1,51))

  %if a channel in a session has more than COVERAGE_THRESHOLD of the time with
  %artifacts, the entire channel gets rejected for that session

  artifact_1_session_cvg_sel = artifact_1_session_cvg(artifact_1_session_cvg.coverage >= ELECTRODE_COVERAGE_THRESHOLD,:);
  artifact_2 = bml_annot_rowbind(artifact_1,artifact_1_session_cvg_sel);
  cfg=[];
  cfg.groupby = 'label';
  artifact_2 = bml_annot_consolidate(cfg,artifact_2);

%   %creating ft_raw from annotations for visualization
%   cfg=[];
%   cfg.template = D_hg_env;
%   cfg.annot_label_colname='label';
%   artifact2_raw = bml_annot2raw(cfg,artifact_2);
% 
%   %raster plot of artifacts for session 1
%   f=figure();
%   bml_plot_raster(artifact2_raw)

  %% cheking coverage per connector group
  %if several channels of the same connector group have an artifact, reject
  %the entire connector group

  %adding connector information to artifac annotation table
  electrode.conn_label = strcat({'conn'},num2str(electrode.connector));
  artifact_2.conn_label = bml_map(artifact_2.label, electrode.electrode, electrode.conn_label);

% 	%calculating absolute value envelope at 1Hz (1s chunks)
%   cfg=[];
%   cfg.freq=ENVELOPE_BIN_SIZE_SECONDS;
%   D_hg_env = bml_envelope_binabs(cfg,D_hg);

  %for each connector and bin, count number of faulty channels
  cfg=[];
  cfg.roi = bml_raw2annot(D_hg_eltype_env);
  cfg.annot_label_colname = 'conn_label';
  connector_artifact_2_cvg_raw = bml_annot2raw(cfg,artifact_2);

%   f=figure();
%   cfg.colorbar = true;
%   bml_plot_raster(cfg,connector_artifact_2_cvg_raw)

  %detecting faulty connectors
  cfg=[];
  cfg.threshold = CONNECTOR_THRESHOLD;
  connector_artifact_3 = bml_annot_detect(cfg,connector_artifact_2_cvg_raw);

  if ~isempty(connector_artifact_3)
    %for each period a connector is faulty, create table with all channels
    %corresponding to that connctor
    cfg=[];
    cfg.groupby_x='conn_label'; %grouping variable in electrode table
    cfg.groupby_y='label'; %corresponding grouping variable in connector_artifact_3
    artifact_4=bml_annot_intersect(cfg,electrode,connector_artifact_3);

    if ~isempty(artifact_4)
      %combining with previously detected artifacts
      artifact_4.label = artifact_4.electrode;
      artifact_5 = bml_annot_rowbind(artifact_2, artifact_4);
      cfg=[];
      cfg.groupby = 'label';
      artifact_5 = bml_annot_consolidate(cfg,artifact_5);
    else
      artifact_5 = artifact_2;
    end
  else
    artifact_5 = artifact_2;
  end

  %final raster plot for artifacts
  cfg=[];
  cfg.template = D_hg_eltype_env;
  cfg.annot_label_colname = 'label';
  artifact_5_raw = bml_annot2raw(cfg,artifact_5);

  clf(f); set(f,'Position',[0 0 600 600]);
  cfg.trial_name='session';
  bml_plot_raster(cfg,artifact_5_raw)
  saveas(f,['figures/' SUBJECT '_' pname '_artifact_mask.png'])

  artifact_5.pname = repmat({pname},height(artifact_5),1);
  
%% saving  artifact annotation table
  artifact = bml_annot_rowbind(artifact,...
    artifact_5(:,{'id','starts','ends','duration','label','conn_label','pname'}));
  
end

cd(PATH_SYNC)
artifact_annot_path = ['annot/' SUBJECT '_artifact_criteria_' CRITERIA '.txt'];


% AM commented out archiving because permissions for sync and annot folders are not available

% % % % %archiving 
% % % % if isfile(artifact_annot_path)
% % % %   copyfile(artifact_annot_path,...
% % % %         [PATH_SYNC filesep 'archive' filesep SUBJECT '_artifact_criteria_' CRITERIA '_' datestr(now,'yyyymmdd_HHMM') '.txt'])
% % % % end





% AM changed write-to folder because permissions for sync and annot folders are not available

bml_annot_write(artifact,['annot/' SUBJECT '_artifact_criteria_' CRITERIA '_denoised.txt']);
bml_annot_write(artifact,[PATH_PROTOCOL, filesep, 'annot', filesep, SUBJECT '_artifact_criteria_' CRITERIA '_denoised.txt']);

