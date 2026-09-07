 %%%% align response timecourses by a sync event, then group together timecourses according to a trial condition

 % inputs:
%       1. trials: table which contain:
%           -EITHER resp_unaligned - a ntrials*1 cell array, with each containing the response of this channel on this trial
%              ....... OR resp_aligned - a ntrials*n_timepoints array... this is the trials.resp_aligned output from align_timecourses
%           -trials.times - must contain the syncing time in each trial.... only required if resp_aligned not provided
%           -a variable with name matching op.time_align_var for syncing responses.... only required if resp_aligned not provided
%       2. op: struct which contains
%            -sort_cond = the name of a variable in trials table, according to which we will group response timecourses
%            -time_align_var = the name of a variable in trials table containing times that responses will be aligned to on each trial.... only required if resp_aligned not provided
%            -sort_cond_vals = list of all values of value of sort_cond to include, in the order you want
%


%       
% outputs: 
%       1. trials_out = original trials table appended with resp_aligned (responses aligned to intratrial event of inerest)
%       2. align_stats = struct with fields containing simple analyses of aligned timecourses, including timecourse mean, sem, sem bar lims (for plotting), timepoints on each side of sync point
%               .... this contains align_stats.times_aligned added - match this with trials.resp_aligned for plotting
%               .... this will be empty if trials.resp_aligned was provided, as align_timecourses.m will be skipped
%       3. resp_grpd = table with row for each value of the op.sort_cond; contains all responses within this condition, as well as mean, std, sem
%       4. cfg_out = original cfg struct plus defaults that were filled in


 function [trials_out, align_stats, resp_grpd, op_out] = sort_responses_by_condition(trials,op)
      field_default('op','sort_cond_vals',{}); 

     if  ~ismember('resp_aligned', trials.Properties.VariableNames)
         [trials, align_stats, op] = align_timecourses(trials, op);
     else
        align_stats = [];
     end

     % process sorting condition info
    assert(isfield(op,'sort_cond') && any(contains(trials.Properties.VariableNames,op.sort_cond)))
    trials.sort_cond = trials{:,op.sort_cond}; 

    if isempty(op.sort_cond_vals) % if not specified, use all sort condition values
        op.sort_cond_vals = unique( trials.sort_cond );
    end


    if isnumeric(op.sort_cond_vals) 
        op.sort_cond_vals = op.sort_cond_vals(~isnan(op.sort_cond_vals)); % remove NaN condition labels
        op.sort_cond_vals = cellstr(string(op.sort_cond_vals)); % turn from numeric to string
    end
    [~,trial_cond_ind] = ismember(string(trials.sort_cond),op.sort_cond_vals); 
    nconds = length(op.sort_cond_vals);

    % make table of sorting conditions
    celcol = cell(nconds,1);
    resp_grpd = table(reshape(op.sort_cond_vals, [], 1), celcol, celcol, 'VariableNames',{'condval','resp','resp_mean'}); 

    for icond = 1:nconds
        these_trial_inds = trial_cond_ind == icond;
        resp_grpd.resp{icond} = trials.resp_aligned(these_trial_inds,:);
        resp_grpd.resp_mean{icond} = mean(resp_grpd.resp{icond},1,'omitnan');
        resp_grpd.std{icond} = std(resp_grpd.resp{icond}, 'omitnan'); % stdev of response timecourses
        resp_grpd.n_good_trials{icond} = sum(~isnan(resp_grpd.resp{icond})); % number of usable trials for this aligned timepoint
        resp_grpd.sem{icond} = resp_grpd.std{icond} ./ sqrt(resp_grpd.n_good_trials{icond});
    end
    
    trials_out = trials; 
    op_out = op; 
 end