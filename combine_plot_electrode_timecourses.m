 %%%% average timecourses of electrodes and plot


% function combine_plot_electrode_timecourses(resp,op)


field_default('op','newfig', true); 
field_default('op','analyze_responsive_elcs_only',1);
field_default('op','analyze_tuned_elcs_only',1);

% if analyze_tuned_elcs_only==true, you must define a parameter op.tuning_param to serve as inclusion criteria
% if cfg.analyze_tuned_elcs_only==true, then only include elcs where op.tuning_param<op.tuning_alpha
field_default('op','analyze_tuned_elcs_only',1);
    field_default('op','tuning_alpha',0.05); 
    
[op,resp] = define_brain_regions(op,resp); 

if op.newfig 
    hfig = figure('color','w','WindowState', 'maximized');
end

% optional - eliminate electrodes marked as nonresponsive
if op.analyze_responsive_elcs_only || ~ismember('rspv',resp.Properties.VariableNames)
    resp = resp(resp.rspv,:); 
end




op.sort_cond = 'learn_con'; 
op.tuning_param = 'p_prod_learn'; 
op.time_align_var = 't_prod_on'; % speech onset




% optional - eliminate electrodes marked as not tuned to parameter of interest
if op.analyze_tuned_elcs_only
    tuned_elcs = resp{:,op.tuning_param} < op.tuning_alpha; 
    resp = resp(tuned_elcs,:); 
end






n_elc = height(resp);
for ielc = 1:n_elc
    subind = find(string(subs.sub) == resp.sub{ielc});
    trials_this_elc = subs.trials{subind}; 
    trials_this_elc.resp_unaligned = resp.timecourse{ielc};
    [trials_this_elc, align_stats, resp_grpd, ~] = sort_responses_by_condition(trials_this_elc,op); 
    
    
    
    %% NEED TO SAVE THE ABOVE OUTPUTS TO A TABLE
    % ..... then sort by region, average and sem over region, and present each region in a subplot
    % ............ can we use plot_resp_timecourse for this?  hopefully can just build in a parameter to skip sort_responses_by_condition if we can provide resp_grpd


end

% end