%% currently only working with Triplet

function P08A09_highgamma_from_denoised(SUBJECT, param)
% protocol P08 in Pitt data, protocol A09 in MGH data


% AM 2022/11/2

% compute high gamma (averaged across 70-150hz) from vibration-denoised trialwise data....
% .......... and save the data



% the vibration-denoised data has already been filtered, so skip filtering step

% % % % % HIGH_PASS_FILTER = 'yes'; %should a high pass filter be applied
% % % % % HIGH_PASS_FILTER_FREQ = 1; %cutoff frequency of high pass filter
% % % % % 
% % % % % do_bsfilter = 'yes'; 
% % % % % line_noise_harm_freqs=[60 120 180 240]; % for notch filters for 60hz harmonics

%% load packages
ft_defaults
bml_defaults
format long

%% Definig paths
%SUBJECT='DBS3022';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_PROTOCOL = 'C:\Users\amsmeier\Documents\MATLAB\P09_artifact_criteria_E';
PATH_SAVE_PREPROCESSED = [PATH_SUBJECT '/Preprocessed Data/FieldTrip'];

cd(PATH_SYNC)

session= bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);

%% Loading FieldTrip data 
load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_raw_filt_trial_denoised.mat'],'D');
nTrials = numel(D.trial);

%remasking nans with zeros
cfg=[];
cfg.value=0;
cfg.remask_nan=true;
D=bml_mask(cfg,D);

%% working in protocol folder
% cd(PATH_PROTOCOL)

%% loading electrode type band table
if ~exist('el_band','var')
  param = readtable('artifact_E_params.txt');
  param_default = param(param.subject == "default",:);
  param_subject = param(strcmp(param.subject,SUBJECT),:);
  if ~isempty(param_subject)
    param = bml_annot_rowbind(param_default(~ismember(param_default.name,param_subject.name),:),param_subject);
  end
end

%%% the vibration-denoised data has already been filtered, so skip filtering step

% % % % % % % % %% Applying High Pass Filter
% % % % % % % % cfg=[];
% % % % % % % % cfg.hpfilter=HIGH_PASS_FILTER;
% % % % % % % % cfg.hpfreq=HIGH_PASS_FILTER_FREQ;
% % % % % % % % cfg.hpfilttype='but';
% % % % % % % % cfg.hpfiltord=5;
% % % % % % % % cfg.hpfiltdir='twopass';
% % % % % % % % cfg.bsfilter=do_bsfilter;
% % % % % % % % cfg.bsfreq= [line_noise_harm_freqs-1; line_noise_harm_freqs+1]'; % notch filters for 60hz harmonics
% % % % % % % % cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'};
% % % % % % % % D_hpf = ft_preprocessing(cfg,D);

%% Artifact rejection
% iterating over bands and electrode types
artifact = table();
f = figure();
for idx = 1:height(param)
  
  fprintf('doing %s %s \n',SUBJECT,param.name{idx});
  
  el_type = strip(param.electrode_type{idx});
    wav_width = param.wav_width(idx);
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
  D_eltype = ft_selectdata(cfg,D);

  if isempty(D_eltype.label)
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


% compute log-spaced frequencies between wav_freq_min and wav_freq_max
    nfreqs = param.n_wav_freqs(idx); 
    nchannels = length(D_eltype.label);
    wav_freqs = round(logspace(log10(param.wav_freq_min(idx)),log10(param.wav_freq_max(idx)),nfreqs));
    D_multifreq_eltype = cell(nfreqs,1);
    
    normed_pow = cell(1,nTrials); 
    med_pow = NaN(nchannels,nTrials,nfreqs);
    for ifreq = 1:nfreqs
      %calculating absolute value envelope at 1Hz (1s chunks)
      cfg=[];
      cfg.out_freq = 100;
      cfg.wav_freq = wav_freqs(ifreq);
      cfg.wav_width = wav_width;
      D_multifreq_eltype{ifreq} = bml_envelope_wavpow(cfg,D_eltype);
      
      
      D_multifreq_eltype{ifreq}.med_pow_per_block = NaN(nchannels, nTrials); % initialize
      for iblock = 1:nTrials % for each block, normalize by median power
        % rows are channels, so take the median across columns (power at timepoints for each channel)
          D_multifreq_eltype{ifreq}.med_pow_per_block(:,iblock) = median(D_multifreq_eltype{ifreq}.trial{iblock},2);
          % normalize power by median values within each channel for this block
          %%% normed_pow will be filled with all normed powers across blocks and frequencies; we will average across the 3rd dimension (frequency)
          normed_pow{iblock}(:,:,ifreq) = D_multifreq_eltype{ifreq}.trial{iblock} ./ D_multifreq_eltype{ifreq}.med_pow_per_block(:,iblock);
      end
      med_pow(:,:,ifreq) = D_multifreq_eltype{ifreq}.med_pow_per_block; % median powers per block/channel in this freq
    end
    
    D_hg_eltype{idx} = struct; % initialize; averaged high gamma
        D_hg_eltype{idx}.hdr = D_multifreq_eltype{1}.hdr;
        D_hg_eltype{idx}.trial = D_multifreq_eltype{1}.trial;
        D_hg_eltype{idx}.trial = cell(1,nTrials); % to be filled
        D_hg_eltype{idx}.time = D_multifreq_eltype{1}.time;
        D_hg_eltype{idx}.label = D_multifreq_eltype{1}.label;
        if isfield(D_multifreq_eltype{1},'sampleinfo')
            D_hg_eltype{idx}.sampleinfo = D_multifreq_eltype{1}.sampleinfo;
        end

    % get averaged high gamma
    med_pow_mean = mean(med_pow,3); % channel/block median powers, averaged across frequencies of interest
    for iblock = 1:nTrials
        D_hg_eltype{idx}.trial{iblock} = mean(normed_pow{iblock},3);
        % multiply the [channel/block]-specific median powers back, to differentiate between absolute power values of channels and trials
        D_hg_eltype{idx}.trial{iblock} = D_hg_eltype{idx}.trial{iblock} .* med_pow_mean(:,iblock); 
    end
    clear D_multifreq_eltype normed_pow
end

non_empty_elc_types = ~cellfun(@isempty,D_hg_eltype); 
cfg = [];
D_hg_trial = ft_appenddata(cfg, D_hg_eltype{non_empty_elc_types});

save([PATH_SAVE_PREPROCESSED filesep SUBJECT '_ft_hg_trial_denoised.mat'],'D_hg_trial','-v7.3');

end
