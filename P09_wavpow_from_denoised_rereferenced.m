%% currently only working with Triplet

% get high gamma timecourses for each trial, using high gamma definition from artifact criterion E
% 
% run this script after p09 (common average rereferencing) has been performed
%
% AM 

function P09_wavpow_from_denoised_rereferenced(op)

% Loading packages
ft_defaults
bml_defaults
format long

%% Definig paths, load annotation tables
setpaths_dbs_triplet()
set_project_specific_variables()
field_default('op','out_freq',100); % downsample rate in hz

% % % Load FieldTrip raw data - artifact-masked and rereferenced
fprintf('\n* Computing %s power (reref method=%s; artifact criteria %s) for subject %s...',op.resp_signal, op.rereference_method, op.art_crit, op.sub)
load([FT_FILE_PREFIX 'raw_filt_trial_ar-',op.art_crit, '_ref-',op.rereference_method, op.denoise_string, '.mat'], 'D_trial_ref');
ntrials_raw = numel(D_trial_ref.trial);

%remasking nans with zeros
cfg1=[];
cfg1.value=0;
cfg1.remask_nan=true;
D_trial_ref=bml_mask(cfg1,D_trial_ref);
 

% get subject-specific artifact criteria
if ~exist('el_band','var') %%%%%%%%%%%% AM 2024/4/20 not sure what this conditional is for - expose contents and remove if statement? 
  param_default = artparam(artparam.subject == "default",:);
  param_subject = artparam(strcmp(artparam.subject,op.sub),:);
  if ~isempty(param_subject)
    artparam = bml_annot_rowbind(param_default(~ismember(param_default.name,param_subject.name),:),param_subject);
  end
end

%% 
% use artifact-detection parameters to re-compute wavpow from rereferenced data
n_eltypes = height(artparam);
for i_eltype = 1:n_eltypes % handle each electrode type 
    cfg1 = [];
    cfg1.out_freq = op.out_freq; 
    cfg1.suppress_output = 1; 
    cfg1.param = artparam(i_eltype,:); 
    D_wavpow_eltype{i_eltype} = multifreq_avg_power(cfg1, D_trial_ref);
end

% combine electrode types
non_empty_elc_types = ~cellfun(@isempty,D_wavpow_eltype); 
cfg1 = [];
D_wavpow_trial = ft_appenddata(cfg1, D_wavpow_eltype{non_empty_elc_types});

save([FT_FILE_PREFIX op.resp_signal '_trial_ar-',op.art_crit, '_ref-',op.rereference_method, op.denoise_string, '.mat'], 'D_wavpow_trial','-v7.3');






