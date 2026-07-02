%% currently only working with Triplet

% get high gamma timecourses for each trial, using high gamma definition from artifact criterion E
% 
% run this script after p09 (common average rereferencing) has been performed
%
% AM 

function P09_highgamma_from_denoised_rereferenced(SUBJECT,ft_file,ARTIFACT_CRIT,savename)

% Loading packages
ft_defaults
bml_defaults
format long

% % % % % Loading parameters
% SUBJECT='DBS3012';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
% ARTIFACT_CRIT = 'E'; 
SAMPLE_RATE = 100; % downsample rate in hz for high gamma traces

%%%%%%%% change this path once crit E is finalized
PATH_ART_PROTOCOL = 'Z:\DBS\Batch\P08_artifact_criteria_E';

cd(PATH_SYNC)

% % % Load annotation tables

electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
% artifact = bml_annot_read(['annot/' SUBJECT '_artifact_criteria_' ARTIFACT_CRIT '.txt']); % use multi-frequency-averaged high gamma
art_param = readtable([PATH_ART_PROTOCOL filesep 'artifact_' ARTIFACT_CRIT '_params.txt']);
empty_electrode = bml_annot_read(['annot/' SUBJECT '_empty_electrode.txt']);
% % % cue_presentation = bml_annot_read(['annot/' SUBJECT '_cue_presentation.txt']);
trials = bml_annot_read(['annot/' SUBJECT '_coding.txt']);


% % % Load FieldTrip raw data - artifact-masked and rereferenced
load(ft_file);
ntrials_raw = numel(D_trial_ref.trial);

%remasking nans with zeros
cfg=[];
cfg.value=0;
cfg.remask_nan=true;
D_trial_ref=bml_mask(cfg,D_trial_ref);
 
%% 
% use artifact-detection parameters to re-compute high gamma
n_eltypes = height(art_param);
for i_eltype = 1:n_eltypes % handle each electrode type 
% compute log-spaced frequencies between wav_freq_min and wav_freq_max

  fprintf('doing %s %s \n',SUBJECT,art_param.name{i_eltype});
  
  el_type = strip(art_param.electrode_type{i_eltype});
  wav_width = art_param.wav_width(i_eltype);
  env_mult_factor =  art_param.env_mult_factor(i_eltype);
  pname = strip(art_param.name{i_eltype});

  % run each electrode type individually, in case different parameters were used during artifact detection
  cfg=[];
  cfg.channel = [el_type,'_*'];
  D_trial_ref_eltype = ft_selectdata(cfg,D_trial_ref); 
  
    nfreqs = art_param.n_wav_freqs(i_eltype); 
    nchannels = length(D_trial_ref_eltype.label);
    wav_freqs = round(logspace(log10(art_param.wav_freq_min(i_eltype)),log10(art_param.wav_freq_max(i_eltype)),nfreqs));
    D_multifreq_eltype = cell(nfreqs,1);
    
    normed_pow = cell(1,ntrials_raw); 
    med_pow = NaN(nchannels,ntrials_raw,nfreqs);
    for ifreq = 1:nfreqs
      %calculating absolute value envelope
      cfg=[];
      cfg.out_freq = SAMPLE_RATE;
      cfg.wav_freq = wav_freqs(ifreq);
      cfg.wav_width = wav_width;
      cmd  = 'D_multifreq_eltype{ifreq} = bml_envelope_wavpow(cfg,D_trial_ref_eltype);';
        evalc(cmd); % use evalc to suppress console output
      
      
      D_multifreq_eltype{ifreq}.med_pow_per_block = NaN(nchannels, ntrials_raw); % initialize
      for iblock = 1:ntrials_raw % for each block, normalize by median power
        % rows are channels, so take the median across columns (power at timepoints for each channel)
          D_multifreq_eltype{ifreq}.med_pow_per_block(:,iblock) = median(D_multifreq_eltype{ifreq}.trial{iblock},2);
          % normalize power by median values within each channel for this block
          %%% normed_pow will be filled with all normed powers across blocks and frequencies; we will average across the 3rd dimension (frequency)
          normed_pow{iblock}(:,:,ifreq) = D_multifreq_eltype{ifreq}.trial{iblock} ./ D_multifreq_eltype{ifreq}.med_pow_per_block(:,iblock);
      end
      med_pow(:,:,ifreq) = D_multifreq_eltype{ifreq}.med_pow_per_block; % median powers per block/channel in this freq
    end
    
    D_hg_eltype{i_eltype} = struct; % initialize; averaged high gamma
        D_hg_eltype{i_eltype}.hdr = D_multifreq_eltype{1}.hdr;
        D_hg_eltype{i_eltype}.trial = D_multifreq_eltype{1}.trial;
        D_hg_eltype{i_eltype}.trial = cell(1,ntrials_raw); % to be filled
        D_hg_eltype{i_eltype}.time = D_multifreq_eltype{1}.time;
        D_hg_eltype{i_eltype}.label = D_multifreq_eltype{1}.label;
   
    % get averaged high gamma
    med_pow_mean = mean(med_pow,3); % channel/block median powers, averaged across frequencies of interest
    for iblock = 1:ntrials_raw
        D_hg_eltype{i_eltype}.trial{iblock} = mean(normed_pow{iblock},3);
        
        % multiply the [channel/block]-specific median powers back, to differentiate between absolute power values of channels
        D_hg_eltype{i_eltype}.trial{iblock} = D_hg_eltype{i_eltype}.trial{iblock} .* med_pow_mean(:,iblock); 
    end
    clear D_multifreq_eltype normed_pow
end

% combine electrode types
cfg = [];
D_hg = ft_appenddata(cfg, D_hg_eltype{:});

save(savename, 'D_hg','-v7.3');






