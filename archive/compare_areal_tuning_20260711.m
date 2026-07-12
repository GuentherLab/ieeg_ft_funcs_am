 %%%% check whether there is a nonrandom distribution of significantly tuned electrodes across areas
 % op.param can either be just the name of the param, or it can be {name, idx}, where idx is the column within the name table variable of the resp table

function subs = compare_areal_tuning(resp,op)

% close all

field_default('op','newfig', true); 
field_default('op','alpha', 0.05); 
field_default('op','bar_face_color', [0.5 0.5 0.5]);
field_default('op','analyze_responsive_elcs_only',1);
field_default('op','warn_about_unassigned_elcs',0); 
field_default('op','param','p_prod'); % this variable should generally be defined by the calling func; project-specific


field_default('op','separate_individual_subs',0); 
    field_default('op','n_subs_per_row',4); 

if any(contains(resp.Properties.VariableNames,'HCPMMP1_label_1')); op.atlas = 'hcp_distal'; end
op = define_brain_regions(op); 


if op.newfig 
    hfig = figure('color','w','WindowState', 'maximized');
end

if op.analyze_responsive_elcs_only || ~ismember('rspv',resp.Properties.VariableNames)
    resp = resp(rows_to_analyze,:); 
end

if iscell(op.param) % if we need to select only 1 column from the table variable
    resp.paramvals = resp{:,op.param{1}}(:,op.param{2});
    field_default('op','full_param_string',[op.param{1}, '_', num2str(op.param{2})]); % display name of the param - default to its name in resp table with index number
else
    resp.paramvals = resp{:,op.param};
    field_default('op','full_param_string',op.param); % display name of the param - default to its name in resp table
end
    
subs = table(sort(unique(resp.sub)),'VariableNames',{'sub'}); 
if op.separate_individual_subs
    nsubs = height(subs);
    ncols = op.n_subs_per_row;
    nrows = ceil([nsubs+1]/ncols);

    for isub = 1:nsubs
        subplot(nrows,ncols,isub)
        thissub = subs.sub{isub}; 
        resp_sub = resp(strcmp(resp.sub,thissub),:); % elcs for this sub
        op.subtitle = thissub; 
        subs.areastats{isub} =  sort_and_plot_elcs(resp_sub,op);
    end

    % plot the grand avg plot across subs
    subplot(nrows,ncols,nsubs+1)
    op.subtitle = 'all_subs'; 
    sort_and_plot_elcs(resp,op);
    
    
elseif ~op.separate_individual_subs
    op.subtitle = 'all_subs'; 
    sort_and_plot_elcs(resp,op)
    areastats = {}; 
end

titlestr = [op.full_param_string];
htitle = sgtitle(titlestr); 

%% 

function areastats = sort_and_plot_elcs(resp_to_plot,op)

    paramvalid = ~isnan(resp_to_plot.paramvals) & resp_to_plot.paramvals ~= 0; % electrodes with usable p values - p must not be zero
    paramsgn = resp_to_plot.paramvals < op.alpha & paramvalid; % analyzable electrodes significantly tuned for param of interest
    
    nelc = height(resp_to_plot);
    resp_to_plot.region = cell(nelc,1); 
    
    nregions = size(op.regiondef,1) + 1;
    areastats = table([op.regiondef(:,1);'all'], [op.regiondef(:,2);{{'all'}}], nan(nregions,2), 'VariableNames', {'region','subareas','ebar_lims'});
    
    natlas = length(op.atlas_var_names); 
    atlastab = table(op.atlas_var_names, cell(natlas,1), nan(natlas,1), nan(natlas,1), 'VariableNames', {'atlas', 'lblcnts', 'n_labeled_elcs', 'n_unlabeled_elcs'});
    
    for iregion = 1:nregions
        thisregion = areastats.region{iregion};

        if thisregion == "all"
            regionmatch = true(nelc,1); 
        else
            regionmatch = false(nelc,1); 
            for iatlas = 1:natlas
                thisatlas = op.atlas_var_names{iatlas};
                if ismember(thisatlas,resp_to_plot.Properties.VariableNames)
                    regionmatch = regionmatch | any( table2array(rowfun(@(x)strcmp(x,areastats.subareas{iregion}),resp_to_plot,'InputVariables',thisatlas)), 2);
                end
            end
        end

        resp_to_plot.region(regionmatch) = repmat({thisregion},nnz(regionmatch),1);    areastats.nelc(iregion) = nnz(regionmatch); % total electrodes in this region
        areastats.nelc_valid(iregion) = nnz(regionmatch & paramvalid); % number of analyzable electrodes in this region for the param of interest 
        areastats.nelc_sgn(iregion) = nnz(regionmatch & paramsgn); % number of analyzable electrodes significantly tuned for param of interest in this region
        areastats.prop_sgn(iregion) = areastats.nelc_sgn(iregion) / areastats.nelc_valid(iregion); % proportion of tuned electrodes in this region
    
        %%%% compute error bar values - 95% confidence intervals using binomial test on each area independently
        binomial_test_op.alpha=.0001:.0001:.9999; 
        p = binocdf(areastats.nelc_sgn(iregion), areastats.nelc_valid(iregion), binomial_test_op.alpha); %%%% note - this should be the matlab inbuilt version of binocdf - not the fieldtrip version - check your path
        areastats.ebar_lims(iregion,1:2) = binomial_test_op.alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]);
    end
    areastats.Properties.RowNames = areastats.region; 
    
    % warn about electrodes that were not assigned a region
    %%%% need to update this so that it only doesn't treat ecog/dbs electrodes as 'unlabeled' for the purposes of dbs/ecog atlases
    all_subareas = cat(2,[areastats.subareas{:}]'); 
    for iatlas = 1:natlas
        thisatlas = op.atlas_var_names{iatlas}; 
        if ismember(thisatlas,resp.Properties.VariableNames)
            [B,BG,BP] = groupcounts(resp{:,thisatlas}); atlastab.lblcnts{iatlas} = table(BG, B, BP./100, 'VariableNames', {'label','count','proportion'}); clear B BG BP
            atlastab.lblcnts{iatlas}.has_region = cellfun(@(x)ismember(x,all_subareas), atlastab.lblcnts{iatlas}.label); 
            missing_reg = ~atlastab.lblcnts{iatlas}.has_region;
            atlastab.n_labeled_elcs(iatlas) = sum(atlastab.lblcnts{iatlas}.count(~missing_reg)); 
            atlastab.n_unlabeled_elcs(iatlas) = sum(atlastab.lblcnts{iatlas}.count(missing_reg)); 
            atlastab.n_elcs(iatlas) = atlastab.n_labeled_elcs(iatlas) + atlastab.n_unlabeled_elcs(iatlas); 
            if op.warn_about_unassigned_elcs && any(missing_reg)
                fprintf([thisatlas ' labels not assigned a region (', num2str(atlastab.n_unlabeled_elcs(iatlas)), '/', num2str(atlastab.n_elcs(iatlas)) '):' ])
                missing_label = atlastab.lblcnts{iatlas}(missing_reg,:)
            end
        end
    end
    
    % get the overall proportion of tuned electrodes with a region assigned and analyzable p values
    %%%%% alternative computation: region_assigned = ~cellfun(@isempty,resp_temp.region); proportion_signficant_overall = mean( paramsgn(paramvalid & region_assigned) ); 
    
    proportion_signficant_overall = sum(areastats.nelc_sgn) / sum(areastats.nelc_valid); 
    
    % number of tuned electrodes per region if they were randomly distributed across regions
    expected_sgn_per_region_random = proportion_signficant_overall * areastats.nelc_valid; 
    
    [chi_significant, chi_p, chi_stats] = chi2gof([1:nregions]', 'Frequency',areastats.nelc_sgn, 'Expected',expected_sgn_per_region_random, 'Emin',0);
    % chi_p
    
    
    %% plotting
    
    hbar = bar(areastats.prop_sgn);
    
    hold on
    
    ebar_neg =  max(areastats.prop_sgn - areastats.ebar_lims(:,1), 0); % proportion should not go below zero 
    ebar_pos =  -areastats.prop_sgn + areastats.ebar_lims(:,2); 
    h_ebar = errorbar([1:nregions]', areastats.prop_sgn, ebar_neg, ebar_pos,'--');
    h_ebar.LineWidth = 0.8;
    h_ebar.LineStyle = 'none';
    h_ebar.Color = [0 0 0];
    
    hax = gca;
    hax.XTickLabels = areastats.region;
    hyline = yline(op.alpha);
    set(0, 'DefaultTextInterpreter', 'none')
    hbar.FaceColor = op.bar_face_color; 

    if isfield(op,'subtitle')
        subtitle([op.subtitle, '... ', num2str(nnz(paramvalid)), ' electrodes'])
    end
    
    box off
    ylabel('fraction of electrodes tuned')
    
    hold off

end
   
end
