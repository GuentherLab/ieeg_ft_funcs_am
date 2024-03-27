 %%%% check whether there is a nonrandom distribution of significantly tuned electrodes across areas
 % before running, need to load resp_temp table and set regiondef and param (project-specfic

% close all

vardefault('alpha', 0.05); 
vardefault('bar_face_color', [0.5 0.5 0.5]);
vardefault('analyze_responsive_elcs_only',1);
vardefault('warn_about_unassigned_elcs',0); 

atlas_var_names = {'HCPMMP1_label_1';'fs_anatomy';'DISTAL_label_1';'MOREL_label_1'}; 

if analyze_responsive_elcs_only || ~ismember('rspv',resp.Properties.VariableNames)
    rows_to_analyze = resp.rspv; 
elseif ~analyze_responsive_elcs_only
    rows_to_analyze = 1:height(resp); 
end
resp_temp = resp(rows_to_analyze,:);
paramvals = paramvals(rows_to_analyze); 

paramvalid = ~isnan(paramvals) & paramvals ~= 0; % electrodes with usable p values
paramsgn = paramvals < alpha & paramvalid; % analyzable electrodes significantly tuned for param of interest

nelc = height(resp_temp);
resp_temp.region = cell(nelc,1); 

nregions = size(regiondef,1);
areastats = table(regiondef(:,1), regiondef(:,2), nan(nregions,2), 'VariableNames', {'region','subareas','ebar_lims'});

natlas = length(atlas_var_names); 
atlastab = table(atlas_var_names, cell(natlas,1), nan(natlas,1), nan(natlas,1), 'VariableNames', {'atlas', 'lblcnts', 'n_labeled_elcs', 'n_unlabeled_elcs'});

for iregion = 1:nregions
    thisregion = areastats.region{iregion};
    regionmatch = false(nelc,1); 
    for iatlas = 1:natlas
        thisatlas = atlas_var_names{iatlas};
        if ismember(thisatlas,resp_temp.Properties.VariableNames)
            regionmatch = regionmatch | any( table2array(rowfun(@(x)strcmp(x,areastats.subareas{iregion}),resp_temp,'InputVariables',thisatlas)), 2);
        end
    end
    resp_temp.region(regionmatch) = repmat({thisregion},nnz(regionmatch),1);    areastats.nelc(iregion) = nnz(regionmatch); % total electrodes in this region
    areastats.nelc_valid(iregion) = nnz(regionmatch & paramvalid); % number of analyzable electrodes in this region for the param of interest 
    areastats.nelc_sgn(iregion) = nnz(regionmatch & paramsgn); % number of analyzable electrodes significantly tuned for param of interest in this region
    areastats.prop_sgn(iregion) = areastats.nelc_sgn(iregion) / areastats.nelc_valid(iregion); % proportion of tuned electrodes in this region

    %%%% compute error bar values - 95% confidence intervals using binomial test on each area independently
    % move fieldtrip version of binocdf to bottom of path so that we use matlab inbuilt version
    oldpath = path; path(oldpath, [PATH_FIELDTRIP_CODE filesep '\external\stats']); clear oldpath 
    binomial_tes_alpha=.0001:.0001:.9999; 
    p = binocdf(areastats.nelc_sgn(iregion), areastats.nelc_valid(iregion), binomial_tes_alpha);
    areastats.ebar_lims(iregion,1:2) = binomial_tes_alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]);
end

% warn about electrodes that were not assigned a region
%%%% need to update this so that it only doesn't treat ecog/dbs electrodes as 'unlabeled' for the purposes of dbs/ecog atlases
all_subareas = cat(2,[areastats.subareas{:}]'); 
for iatlas = 1:natlas
    thisatlas = atlas_var_names{iatlas}; 
    if ismember(thisatlas,resp_temp.Properties.VariableNames)
        [B,BG,BP] = groupcounts(resp_temp{:,thisatlas}); atlastab.lblcnts{iatlas} = table(BG, B, BP./100, 'VariableNames', {'label','count','proportion'}); clear B BG BP
        atlastab.lblcnts{iatlas}.has_region = cellfun(@(x)ismember(x,all_subareas), atlastab.lblcnts{iatlas}.label); 
        missing_reg = ~atlastab.lblcnts{iatlas}.has_region;
        atlastab.n_labeled_elcs(iatlas) = sum(atlastab.lblcnts{iatlas}.count(~missing_reg)); 
        atlastab.n_unlabeled_elcs(iatlas) = sum(atlastab.lblcnts{iatlas}.count(missing_reg)); 
        atlastab.n_elcs(iatlas) = atlastab.n_labeled_elcs(iatlas) + atlastab.n_unlabeled_elcs(iatlas); 
        if warn_about_unassigned_elcs && any(missing_reg)
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
if show_barplot

    if newfig 
        hfig = figure('color','w');
    end
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
    hyline = yline(alpha);
    set(0, 'DefaultTextInterpreter', 'none')
    titlestr = [full_param_string, '..... p = ' num2str(chi_p)] ;
%     htitle = title(titlestr); 

    hbar.FaceColor = bar_face_color; 

end

box off
ylabel('fraction of electrodes tuned')

hold off
