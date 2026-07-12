 %%%% average timecourses of electrodes and plot


function combine_plot_electrode_timecourses(resp,op)


field_default('op','newfig', true); 
field_default('op','analyze_responsive_elcs_only',1);

% if analyze_tuned_elcs_only==true, you must define a parameter op.tuning_param to serve as inclusion criteria
% if cfg.analyze_tuned_elcs_only==true, then only include elcs where op.tuning_param<op.tuning_alpha
field_default('op','analyze_tuned_elcs_only',1);
    field_default('op','tuning_alpha'); 
    
[op,resp] = define_brain_regions(op,resp); 

if op.newfig 
    hfig = figure('color','w','WindowState', 'maximized');
end

% optional - eliminate electrodes marked as nonresponsive
if op.analyze_responsive_elcs_only || ~ismember('rspv',resp.Properties.VariableNames)
    resp = resp(resp.rspv,:); 
end

% optional - eliminate electrodes marked as not tuned to parameter of interest
if op.analyze_tuned_elcs_only
    tuned_elcs = resp{op.tuning_param,:} < op.tuning_alpha; 
    resp = resp(tuned_elcs,:); 
end

