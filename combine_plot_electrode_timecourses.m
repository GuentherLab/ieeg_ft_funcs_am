 %%%% average timecourses of electrodes and plot

 %%%%%%%% inputs %%%%%%%%
 % 1. resp - table of electrode responses containing:
 %      EITHER timecourse OR timecourses_aligned
 %      sub
 %      chan
 %      rspv [required if analyze_responsive_elcs_only==1]
 %      [tuning param - required if analyze_tuned_elcs_only==1]
 %      [anatomical area label.... probably HCPMMP1_label_1 or DISTAL_label_1      
 %
 % 2. op [optional]
 %
 %
 %
 %       Recommended: the resp table may montain the variable 'timecourses_aligned', in which alignment will have already been done.
%       This table must have a rowname matching op.time_align_var ....
%            ..... and variable names 'trials' and 'align_stats', which were the outputs of align_timecourses.m
%        If this table is present, then align_timecourses will be skipped here (potentially saves a lot of time)


 function [cond_elc_rgn, align_stats_rgn, resp_grpd_rgn, cfg_rgn] = combine_plot_electrode_timecourses(resp,subs,op)

vardefault('op',struct);
field_default('op','newfig', true); 
field_default('op','analyze_responsive_elcs_only',1);
field_default('op','analyze_tuned_elcs_only',1);
field_default('op','row_col_ratio',4/3); % approx ratio of rows to columns
field_default('op','xline_events',{}); 
field_default('op','sort_cond_vals',{}); 

% if analyze_tuned_elcs_only==true, you must define a parameter op.tuning_param to serve as inclusion criteria
% if cfg.analyze_tuned_elcs_only==true, then only include elcs where op.tuning_param<op.tuning_alpha
field_default('op','analyze_tuned_elcs_only',1);
    field_default('op','tuning_alpha',0.05); 
   



% optional - eliminate electrodes marked as nonresponsive
if op.analyze_responsive_elcs_only || ~ismember('rspv',resp.Properties.VariableNames)
    resp = resp(resp.rspv,:); 
end

% optional - eliminate electrodes marked as not tuned to parameter of interest
if op.analyze_tuned_elcs_only
    tuned_elcs = resp{:,op.tuning_param} < op.tuning_alpha; 
    resp = resp(tuned_elcs,:); 
end

[op,resp] = define_brain_regions(op,resp); 

% make a copy of electrodes table for aligning and averaging responses
n_elc = height(resp);
resp_align = resp(:,{'sub','chan','region'}); 
resp_align.times_aligned = cell(n_elc,1); 
resp_align.resp_grpd = cell(n_elc,1); 

temptab = table; 
for ielc = 1:n_elc % % this loop can take a minute for >1000 electrodes if resp.timecourses_unaligned was not provided
    subind = find(string(subs.sub) == resp.sub{ielc});
    trials_this_elc = subs.trials{subind}; 

    if ismember('timecourses_aligned', resp.Properties.VariableNames)
        trials_this_elc.resp_aligned = resp.timecourses_aligned{ielc}{op.time_align_var,'trials'}{1}.resp_aligned; 
    elseif ~ismember('timecourses_aligned', resp.Properties.VariableNames)
        trials_this_elc.resp_unaligned = resp.timecourse{ielc};
    end


    [~, align_stats_this_elc, resp_grpd, ~] = sort_responses_by_condition(trials_this_elc,op); 
    resp_align.times{ielc} = align_stats_this_elc.times_aligned; 
    resp_align.resp_grpd{ielc} = resp_grpd; 

    temptab = [temptab; repelem(resp_align(ielc, {'sub','chan','region','times'}), height(resp_grpd), 1)];
    
end

% create a table with one row per condition per electrode, with the (within-electrode) aligned mean response on a single condition
%%% we will treat each electrode's mean response as a 'trial' for the purposes of the sort_responses_by_condition function
cond_elc_resp_align = subsref(vertcat(resp_align.resp_grpd{:}), struct('type', '()', 'subs', {{':', {'condval', 'resp_mean','n_good_trials'}}}));
cond_elc_resp_align = renamevars(cond_elc_resp_align, {'resp_mean','condval'}, {'resp_unaligned',op.sort_cond}); % no longer considered 'aligned' because we're now aligning across elcs
cond_elc_resp_align = [temptab, cond_elc_resp_align];
cond_elc_resp_align{:,op.time_align_var} = zeros(height(cond_elc_resp_align),1); % times are already zeroed within trial, so align times are all zero



% align responses to each condition within each region
region_resp = op.regiondef; 

%% plot averaged timecourse within each region

if op.newfig 
    hfig = figure('color','w','WindowState', 'maximized');
end

if op.analyze_tuned_elcs_only
    figtitle = {['electrodes meeting (', op.tuning_param, ' < ', num2str(op.tuning_alpha),')'],''};
else
    figtitle = {'no electrode tuning criteria',''};
end

sgtitle(figtitle);

% % get number of rows and columns
r = 1:op.nregions;  c = ceil(op.nregions./r);
[~, idx] = min(abs(c./r - 1/op.row_col_ratio)); 
op.n_plot_rows = r(idx); op.n_plot_cols = c(idx);


%% construct trial table of xline event variables across subjects
% for the purposes of plotting xline events
subs_plotted = subs(contains(subs.sub, unique(resp.sub)),:);
trials_all_subs = table; 
for isub = 1:height(subs_plotted)
    trials_all_subs = [trials_all_subs; subs_plotted.trials{isub}(:,{op.sort_cond,op.time_align_var,op.xline_events{1,:}})];
end

%% do plotting
cond_elc_rgn = cell(op.nregions,1);
align_stats_rgn = cell(op.nregions,1);
resp_grpd_rgn  = cell(op.nregions,1);
cfg_rgn   = cell(op.nregions,1);
for iregion = 1:op.nregions
    thisregion = region_resp.region{iregion}; 
    hsubplot(iregion) = subplot(op.n_plot_rows,op.n_plot_cols,iregion);

    if region_resp.n_elcs(iregion) > 1 % skip the region if there's less than 2 electrodes to plot
        cond_elc_rgn{iregion} = cond_elc_resp_align(strcmp(cond_elc_resp_align.region,thisregion), :); % make table with only elcs in this region
    
        cfg = [];
        cfg.sort_cond = op.sort_cond; 
        cfg.time_align_var = op.time_align_var;
        cfg.sort_cond_vals = op.sort_cond_vals;
        [cond_elc_rgn{iregion}, align_stats_rgn{iregion}, resp_grpd_rgn{iregion}, cfg_rgn{iregion}] = sort_responses_by_condition(cond_elc_rgn{iregion},cfg);
    
        cfg = [];
        cfg = op;
        cfg.samp_period = cfg_rgn{iregion}.samp_period; 
        cfg.do_condition_sorting = 0; % skip sorting, it's already done here
        cfg.align_stats = align_stats_rgn{iregion}; 
        cfg.resp_grpd = resp_grpd_rgn{iregion}; 
        cfg.newfig = 0; 
        plot_resp_timecourse(trials_all_subs,cfg); % nb: trials table here is only used to plot xline, not to compute response timecourses
    end

    title([thisregion, ' (n=',num2str(region_resp.n_elcs(iregion)),' elc)'])

end



end