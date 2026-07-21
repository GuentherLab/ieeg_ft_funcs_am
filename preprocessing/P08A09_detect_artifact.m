% Detect artifacts, create an artifact annotation table

function [D_masked, thresh_artifact_table] = P08A09_detect_artifact(D,op)

%% load packages
% ft_defaults
% bml_defaults
format long

%% Defining paths, loading artifact parameters
vardefault('op',struct); % initialize options if not present
% field_default('op','sub','DM1007')
% field_default('op','art_crit','E'); % 'E' = 70-250hz high gamma; 'F' = beta; 'G' = other Rohan criterion? 
field_default('op','out_freq',100); % downsample rate in hz


%% Loading FieldTrip data 
ntrials = numel(D.trial);

%remasking nans with zeros
% cfg=[];
% cfg.value=0;
% cfg.remask_nan=true;
% D=bml_mask(cfg,D);


%% loading electrode type band table
% % % % % % % % % % % % % % % % artparam_default = artparam(artparam.subject == "default",:);
% % % % % % % % % % % % % % % % artparam_subject = artparam(strcmp(artparam.subject,op.sub),:);
% % % % % % % % % % % % % % % % if ~isempty(artparam_subject)
% % % % % % % % % % % % % % % %     artparam = bml_annot_rowbind(artparam_default(~ismember(artparam_default.name,artparam_subject.name),:),artparam_subject);
% % % % % % % % % % % % % % % % end


%% Artifact rejection - ECoG channels 

artparam = op.artparam;
electrodes = op.electrodes; 

% iterating over bands and electrode types
thresh_artifact_table = table();
hfig = figure();
for idx = 1:height(artparam)
  
%   fprintf('doing %s %s \n',op.sub,artparam.name{idx});
  pname = strip(artparam.name{idx});
    
    cfg         = [];
    cfg.channel = [artparam.electrode_type{idx}, '*'];
    D_eltype = ft_selectdata(cfg, D);

    if isempty(D_eltype)
        %channel type not available
        continue
    end

  ENVELOPE_BIN_SIZE_SECONDS = artparam.env_bin_size(idx); %envelope bin size in seconds
  THRESHOLD_STD_FACTORS = [artparam.th_factor_std_low(idx), artparam.th_factor_std_high(idx)]; %factors to determine detection thresholds 
  THRESHOLD_FIX = [artparam.th_fix_min(idx), artparam.th_fix_max(idx)]; %fix thresholds to filter data before applying robust estimates
  CONSOLIDATION_TIME_TOLERANCE = artparam.th_consolidation(idx); %min time allowed between consecutive artifacts
  ELECTRODE_COVERAGE_THRESHOLD = artparam.th_frac_coverage(idx); %max allowed fraction of time with artifacts
  CONNECTOR_THRESHOLD = [artparam.th_conn_low(idx), artparam.th_conn_high(idx)]; %detection threshold for number of electrodes in a connector  

  cfg=[];
  cfg.freq = 1/ENVELOPE_BIN_SIZE_SECONDS;
  D_eltype_env = bml_envelope_binabs(cfg,D_eltype);

  %calculating log10 transform (envelopes have log normal distributions)
  D_eltype_env_log10 = bml_apply(@(x) artparam.env_mult_factor(idx) .* log10(x),D_eltype_env);
  
  cfg=[];
  cfg.remask_inf=true;
  cfg.value=NaN;
  D_eltype_env_log10 = bml_mask(cfg,D_eltype_env_log10);
  
  %calculating distribution robust statistics. 
  THRESHOLD = nan(ntrials,2);
  max_v=nan(1,ntrials);
  min_v=nan(1,ntrials);
  for itrial=1:ntrials
    v = reshape(D_eltype_env_log10.trial{itrial},1,[]);
    v1 = v((v>THRESHOLD_FIX(1)) & (v<THRESHOLD_FIX(2)));
    m = median(v1);
    std = bml_robust_std(v1);
    if ~isempty(v1)
      max_v(itrial)=max(v);
      min_v(itrial)=min(v);
      THRESHOLD(itrial,:) = m + THRESHOLD_STD_FACTORS.*std;
    end
  end

  %plotting histogram to assess threshold levels
  if ~isempty(v1)
      clf(hfig); set(hfig,'Position',[0 0 600 600]);
      for itrial=1:ntrials
        subplot(ceil(ntrials/2),2,itrial)
        hold on;
        h=histogram(D_eltype_env_log10.trial{itrial},linspace(min(min_v),max(max_v),61),...
          'FaceAlpha',0.1,'EdgeAlpha',1);
        maxBinCount = max(h.BinCounts);
        plot([THRESHOLD(itrial,1),THRESHOLD(itrial,1)],[0,maxBinCount .* 1.1]);
        plot([THRESHOLD(itrial,2),THRESHOLD(itrial,2)],[0,maxBinCount .* 1.1]);
        %set(gca,'YScale','log')
        title(['session ' num2str(itrial)]);
      end
%       saveas(hfig,[PATH_FIGURES filesep op.sub '_' pname '_artifact_env_log10_hist.png'])
  elseif isempty(v1)
      warning(['For electrode type ''', artparam.electrode_type{idx}, ''' (sub ', op.sub, '), no timepoints found between low threshold (', ...
          num2str(THRESHOLD_FIX(1)), ') and high threshold (', num2str(THRESHOLD_FIX(2)), ')'])
  end

  %detecting segments of time for each channel above threshold
  artifact_eltype_1 = table();
  for itrial=1:ntrials
    cfg=[];
    cfg.threshold = THRESHOLD(itrial,:);
    cfg.trials = itrial;
    artifact_eltype_1 = bml_annot_rowbind(artifact_eltype_1, bml_annot_detect(cfg,D_eltype_env_log10));
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
  [D_eltype_sample, epoch_avgpow_eltype_sample] = bml_redefinetrial(cfg,D_eltype);
 
  D_p = D_eltype_sample;
  E_p = epoch_avgpow_eltype_sample;
  nx=10; ny=floor(numel(D_p.trial)/nx);
  if ny==0
    ny=1; nx=numel(D_p.trial);
  end
  clf(hfig); set(hfig,'Position',[0 0 nx*200 ny*200]);
  for itrial=1:ny
      for j=1:nx
          pidx = (itrial-1)*nx+j;
          l = E_p.label(pidx);
          l_idx = bml_getidx(l,D_p.label);
          subplot(ny,nx,pidx);
          plot(D_p.time{pidx},D_p.trial{pidx}(l_idx,:));
          title(E_p.label(pidx));
      end
  end
%   saveas(hfig,[PATH_FIGURES filesep op.sub '_' pname '_artifact_crit-', op.art_crit,'_snippets_not_denoised.png'])

  
  %consolidating annotations with CONSOLIDATION_TIME_TOLERANCE margin of overlap
  cfg=[];
  cfg.criterion = @(x) (x.starts(end) - max(x.ends(1:(end-1))) < CONSOLIDATION_TIME_TOLERANCE);
  cfg.groupby = 'label';
  artifact_eltype_2 = bml_annot_consolidate(cfg,artifact_eltype_1);

%   %creating ft_raw from annotations for visualization
%   cfg=[];
%   cfg.template = D_hpf_eltype_env_log10;
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
%   D_hpf_eltype_mask = bml_mask(cfg, D_hpf_eltype);
% 
%   cfg=[];
%   cfg.viewmode = 'vertical';
%   cfg.blocksize = 30;
%   cfg.ylim = 'maxmin';
%   cfg.continuous = 'yes';
%   ft_databrowser(cfg,D_hpf_eltype_mask);

  %% rejecting faulty channels 

  %decide which artifacts to include. Usually just ECoG artifacts
  %artifact = bml_annot_rowbind(artifact_ecog_3,artifact_macro_3,artifact_dbs_3);
  artifact_1 = artifact_eltype_2;

  

  cfg = [];
  run_dur_table = table(1, D.time{1}(1), D.time{end}(end), [D.time{end}(end), D.time{1}(1)], 'VariableNames',...
                      {'id', 'starts',      'ends',             'duration'}); 
  cfg.groupby = 'label';
  artifact_1_session_cvg = bml_annot_coverage(cfg,artifact_1,run_dur_table);

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
%   cfg.template = D_hpf_env;
%   cfg.annot_label_colname='label';
%   artifact2_raw = bml_annot2raw(cfg,artifact_2);
% 
%   %raster plot of artifacts for session 1
%   f=figure();
%   bml_plot_raster(artifact2_raw)

  %% checking coverage per connector group
  %if several channels of the same connector group have an artifact, reject
  %the entire connector group

  %adding connector information to artifact annotation table
  electrodes.conn_label = strcat({'conn'},num2str(electrodes.connector));
  artifact_2.conn_label = bml_map(artifact_2.label, electrodes.name, electrodes.conn_label);

% 	%calculating absolute value envelope at 1Hz (1s chunks)
%   cfg=[];
%   cfg.freq=ENVELOPE_BIN_SIZE_SECONDS;
%   D_hpf_env = bml_envelope_binabs(cfg,D_hpf);

  %for each connector and bin, count number of faulty channels
  cfg=[];
  cfg.roi = bml_raw2annot(D_eltype_env);
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
    artifact_4=bml_annot_intersect(cfg,electrodes,connector_artifact_3);

    if ~isempty(artifact_4)
      %combining with previously detected artifacts
      artifact_4.label = artifact_4.name;
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
  cfg.template = D_eltype_env;
  cfg.annot_label_colname = 'label';
  artifact_5_raw = bml_annot2raw(cfg,artifact_5);

%   clf(hfig); 
figure
  set(hfig,'Position',[0 0 600 600]);
  cfg.trial_name='session';
  bml_plot_raster(cfg,artifact_5_raw)
%   saveas(hfig,[PATH_FIGURES filesep op.sub '_' pname '_artifact_mask_not_denoised.png'])

  artifact_5.pname = repmat({pname},height(artifact_5),1);
  
%% consolidate  artifact annotation table
  thresh_artifact_table = bml_annot_rowbind(thresh_artifact_table,...
    artifact_5(:,{'id','starts','ends','duration','label','conn_label','pname'}));
  
end

% apply the mask; set artifactual windows to nan
    cfg2 = [];
    cfg2.annot = thresh_artifact_table; % set manual artifacts to nan
    cfg2.complete_trial = false; % full run data, so can't mask individual trials
    cfg2.value = NaN; 
    cfg2.label_colname = 'label';
    D_masked = bml_mask(cfg2, D); 

% % % % % % % % if an artifact table was already created for this subject with this filename, archive the old file before saving this new file
% % % % % % % if isfile(ARTIFACT_FILENAME_SUB)
% % % % % % %     mkdir([PATH_ANNOT filesep 'archive'])
% % % % % % %   copyfile(ARTIFACT_FILENAME_SUB,...
% % % % % % %         [PATH_ANNOT filesep 'archive' filesep op.sub '_artifact-criteria-' op.art_crit '_not-denoised_' datestr(now,'yyyymmdd_HHMM') '.tsv'])
% % % % % % % end




