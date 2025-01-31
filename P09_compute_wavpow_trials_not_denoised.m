
% get wavelet-transformed power timecourses for each trial, using freq bands definitions from artifact criterion E
% 
% run this script after artifact detection and common average rereferencing has been performed
%
% AM 

function P09_compute_wavpow_trials_not_denoised(op)

% Loading packages
ft_defaults
bml_defaults
format long

%% Defining paths, loading artifact parameters
field_default('op','out_freq',100); % downsample rate in hz

set_project_specific_variables(); % set paths etc. based on data collection site

fieldtrip_savename = [FT_FILE_PREFIX, op.resp_signal, '-trial_ar-',op.art_crit, '_ref-',op.rereference_method,  '_not-denoised.mat'];

% % % Load FieldTrip raw data - artifact-masked and rereferenced
load([FT_FILE_PREFIX 'raw-filt-trial_ar-',op.art_crit, '_ref-',op.rereference_method,  '_not-denoised.mat']);
ntrials_raw = numel(D_trial_ref.trial);

%remasking nans with zeros
cfg=[];
cfg.value=0;
cfg.remask_nan=true;
D_trial_ref=bml_mask(cfg,D_trial_ref);
 
%% 
% use artifact-detection parameters to re-compute high gamma
n_eltypes = height(artparam);
D_avgpow_eltype = cell(n_eltypes,1); 
for i_eltype = 1:n_eltypes % handle each electrode type  
  fprintf('doing %s %s \n',op.sub,artparam.name{i_eltype});

    cfg = [];
    cfg.out_freq = op.out_freq; 
    cfg.suppress_output = 1; 
    cfg.param = artparam(i_eltype,:); 
    D_avgpow_eltype{i_eltype} = multifreq_avg_power(cfg, D_trial_ref); % compute log-spaced frequencies between wav_freq_min and wav_freq_max
    
    if isempty(D_avgpow_eltype{i_eltype})
        %channel type not available
        continue
    end
end

% combine electrode types
cfg = [];
D_avgpow_eltype = D_avgpow_eltype(~cellfun(@isempty,D_avgpow_eltype)); % delete empty rows
D_wavpow = ft_appenddata(cfg, D_avgpow_eltype{:});

save(fieldtrip_savename, 'D_wavpow', '-v7.3');






