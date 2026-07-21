 % ALIGN_TIMECOURSES - realign responses to that they are synced on a specific event (e.g. speech onset) in each trial
%     [trials, align_stats, cfg_out] = align_timecourses(trials, cfg)
%
% inputs:
%       1. trials: table which must contain:
%           -trials.resp_unaligned - a ntrials*1 cell array, with each containing the response of this channel on this trial
%           -trials.times - must contain the syncing time in each trial
%           -a variable with name matching op.time_align_var for syncing responses
%       2. op: struct which must contain
%            time_align_var, which is the name of a variable in trials table containing times that responses will be aligned to on each trial
%
%       
% outputs: 
%       1. trials_out = original trials table appended with resp_aligned (responses aligned to intratrial event of inerest)
%       2. align_stats = struct with fields containing simple analyses of aligned timecourses, including timecourse mean, sem, sem bar lims (for plotting), timepoints on each side of sync point
%               .... this contains align_stats.times_aligned added - match this with trials.resp_aligned for plotting
%       3. cfg_out = original cfg struct plus defaults that were filled in


function [trials_out, align_stats, cfg_out] = align_timecourses(trials, cfg)

% this option specifies how we decide the duration of aligned trials, which required cutting off some long trials and padding shorting trials
field_default('cfg','trial_time_adj_method','median_plus_sd'); % options: mean, median, median_plus_stdev

assert(isfield(cfg,'time_align_var') && any(contains(trials.Properties.VariableNames,cfg.time_align_var)))

% remove trials with no response data
non_empty_trials = ~cellfun(@isempty, trials.resp_unaligned);
trials = trials(non_empty_trials,:);    

trials.align_time = trials{:,cfg.time_align_var}; 

ntrials = height(trials);
 nans_tr = nan(ntrials,1); 

 trials = [trials, table(nans_tr,              nans_tr,...
     'VariableNames',     {'tpoints_pre_onset', 'tpoints_post_onset'})];
 
 % compute sampling interval... assumes a static sample rate
 cfg.samp_period = 1e-5 * round(1e5 * diff(trials.times{1}(1:2))); 
 
 %%% find trial lengths pre- and post- the alignment time
for itrial = 1:ntrials
    % n timepoints before or at align_time
    trials.tpoints_pre_onset(itrial) = nnz(trials.times{itrial} <= trials.align_time(itrial,1)); 
    % n timepoints after align_time
    trials.tpoints_post_onset(itrial) = nnz(trials.times{itrial} > trials.align_time(itrial,1)); 
end
 
% pad or cut each trial to fit a specific size, so that we can align and average trials
align_stats = struct; 
switch cfg.trial_time_adj_method
    case 'median_plus_sd'
        align_stats.n_tpoints_pre_fixed = round(median(trials.tpoints_pre_onset) + std(trials.tpoints_pre_onset)); 
        align_stats.n_tpoints_post_fixed = round(median(trials.tpoints_post_onset) + std(trials.tpoints_post_onset)); 
    case 'median'
        align_stats.n_tpoints_pre_fixed = median(trials.tpoints_pre_onset); 
        align_stats.n_tpoints_post_fixed = median(trials.tpoints_post_onset); 
    case 'max'
        align_stats.n_tpoints_pre_fixed = max(trials.tpoints_pre_onset); 
        align_stats.n_tpoints_post_fixed = max(trials.tpoints_post_onset); 
end
tpoints_tot = align_stats.n_tpoints_pre_fixed + align_stats.n_tpoints_post_fixed; 

% nan-pad or cut trial windows so that they are all the same duration
%%% pad and cut values must be non-negative

resp_aligned = NaN(ntrials, tpoints_tot); % aligned responses for this electrode; rows = trials, columns = timepoints
for itrial = 1:ntrials
   pre_pad = max([0, align_stats.n_tpoints_pre_fixed - trials.tpoints_pre_onset(itrial)]); 
   pre_cut = max([0, -align_stats.n_tpoints_pre_fixed + trials.tpoints_pre_onset(itrial)]); 
   pre_inds = 1+pre_cut:trials.tpoints_pre_onset(itrial); % inds from trials.resp_unaligned... if pre_cut > 0, some timepoints from this trial will not be used
   resp_aligned(itrial, pre_pad+1 : align_stats.n_tpoints_pre_fixed) = trials.resp_unaligned{itrial}(pre_inds); % fill in pre-onset data... fill in electrode responses starting after the padding epoch

   post_pad = max([0, align_stats.n_tpoints_post_fixed - trials.tpoints_post_onset(itrial)]);
   post_cut = max([0, -align_stats.n_tpoints_post_fixed + trials.tpoints_post_onset(itrial)]); 
   post_inds = trials.tpoints_pre_onset(itrial) +  [1 : trials.tpoints_post_onset(itrial)-post_cut]; % inds from trials.resp_unaligned
   resp_aligned(itrial, align_stats.n_tpoints_pre_fixed+1:end-post_pad) = trials.resp_unaligned{itrial}(post_inds); % fill in post-onset data
   
   trials.trial_onset_adjust(itrial) = cfg.samp_period * [pre_pad - pre_cut]; % number of timepoints to add to time landmarks
end
align_stats.mean = mean(resp_aligned,'omitnan'); % mean response timecourse
align_stats.std = std(resp_aligned, 'omitnan'); % stdev of response timecourses
align_stats.std_lims = [align_stats.mean + align_stats.std; align_stats.mean - align_stats.std]; 
align_stats.n_nonnan_trials = sum(~isnan(resp_aligned)); % number of usable trials for this aligned timepoint
align_stats.sem = align_stats.std ./ sqrt(align_stats.n_nonnan_trials);
align_stats.sem_lims = [align_stats.mean + align_stats.sem; align_stats.mean - align_stats.sem]; 

trials.resp_aligned = resp_aligned; 

align_stats.times_aligned = 0.5 + [linspace(-align_stats.n_tpoints_pre_fixed, -1, align_stats.n_tpoints_pre_fixed), linspace(0, align_stats.n_tpoints_post_fixed-1, align_stats.n_tpoints_post_fixed)];
align_stats.times_aligned = cfg.samp_period * align_stats.times_aligned; 

trials_out = trials; 
cfg_out = cfg; 

