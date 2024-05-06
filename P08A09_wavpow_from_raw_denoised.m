%% currently only working with Triplet... need move certain parameters into the calling script to make more generic
% protocol P08 in Pitt data, protocol A09 in MGH data
% compute power in freq band of interest from vibration-denoised trialwise data and save the data



function P08A09_wavpow_from_raw_denoised(op)

%% load packages
ft_defaults
bml_defaults
format long


%% Definig paths
setpaths_dbs_triplet()
set_project_specific_variables()
field_default('op','out_freq',100); % downsample rate in hz

%% Loading FieldTrip data 
load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep op.sub '_ft_raw_filt_trial_denoised.mat'],'D');

%remasking nans with zeros
cfg1=[];
cfg1.value=0;
cfg1.remask_nan=true;
D=bml_mask(cfg1,D);

% get subject-specific artifact criteria
if ~exist('el_band','var') %%%%%%%%%%%% AM 2024/4/20 not sure what this conditional is for - expose contents and remove if statement? 
  param_default = artparam(artparam.subject == "default",:);
  param_subject = artparam(strcmp(artparam.subject,op.sub),:);
  if ~isempty(param_subject)
    artparam = bml_annot_rowbind(param_default(~ismember(param_default.name,param_subject.name),:),param_subject);
  end
end

%%% the vibration-denoised data has already been filtered, so skip filtering step

% iterating over bands and electrode types

% % % % % % % % % % % % % % % f = figure();
n_eltypes = height(artparam); 
for i_eltype = 1:n_eltypes
    cfg1 = [];
    cfg1.out_freq = op.out_freq; 
    cfg1.suppress_output = 1; 
    cfg1.param = artparam(i_eltype,:); 
    D_wavpow_eltype{i_eltype} = multifreq_avg_power(cfg1, D);
end

non_empty_elc_types = ~cellfun(@isempty,D_wavpow_eltype); 
cfg1 = [];
D_wavpow_trial = ft_appenddata(cfg1, D_wavpow_eltype{non_empty_elc_types});

save([FT_FILE_PREFIX, op.resp_signal, '_trial', op.denoise_string, '.mat'],'D_wavpow_trial','-v7.3');

end
