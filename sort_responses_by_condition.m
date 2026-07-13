 %%%% align response timecourses by a sync event, then group together timecourses according to a trial condition

 % inputs:
%       1. trials: table which must contain:
%           -resp_unaligned - a ntrials*1 cell array, with each containing the response of this channel on this trial
%           -a variable with name matching op.time_align_var for syncing responses
%       2. op: struct which must contain
%            -time_align_var = the name of a variable in trials table containing times that responses will be aligned to on each trial
%            -sort_cond = the name of a variable in trials table, according to which we will group response timecourses
%
%       
% outputs: 
%       1. trials_out = original trials table appended with resp_aligned (responses aligned to intratrial event of inerest)
%       2. align_stats = struct with fields containing simple analyses of aligned timecourses, including timecourse mean, sem, sem bar lims (for plotting), timepoints on each side of sync point
%               .... this contains align_stats.xtime added - match this with trials.resp_aligned for plotting
%       3. resp_grpd = table with row for each value of the op.sort_cond; contains all responses within this condition, as well as mean, std, sem
%       4. cfg_out = original cfg struct plus defaults that were filled in


 function [trials_out, align_stats, resp_grpd, op_out] = sort_responses_by_condition(trials,op)

     [trials, align_stats, op] = align_timecourses(trials, op); 


    assert(isfield(op,'sort_cond') && any(contains(trials.Properties.VariableNames,op.sort_cond)))
    trials.sort_cond = trials{:,op.sort_cond}; 


    [unq_conds, ~, trial_cond_ind] = unique( trials.sort_cond );
    if isnumeric(unq_conds) % remove NaN condition labels
        [unq_conds, ~, trial_cond_ind] = unq_conds(~isnan(unq_conds));
        unq_conds = cellstr(num2str(unq_conds));
    end
    nconds = length(unq_conds);

    celcol = cell(nconds,1);
    resp_grpd = table(unq_conds,celcol,celcol,'VariableNames',{'condval','resp','resp_mean'}); 

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