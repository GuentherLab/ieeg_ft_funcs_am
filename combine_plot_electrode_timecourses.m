function [cond_elc_rgn, align_stats_rgn, resp_grpd_rgn, cfg_rgn] = ...
    combine_plot_electrode_timecourses(resp, subs, op)
    
    %% COMBINE_PLOT_ELECTRODE_TIMECOURSES
    % Main function for visualizing pre-warped electrode timecourses
    % Assumes data is already warped: resp.timecourse contains ntrial x ntimepoint responses
    
    %% Set defaults
    vardefault('op', struct)
    field_default('op', 'newfig', true)
    field_default('op', 'analyze_responsive_elcs_only', 1)
    field_default('op', 'analyze_tuned_elcs_only', 1)
    field_default('op', 'tuning_alpha', 0.05)
    field_default('op', 'row_col_ratio', 4/3)
    field_default('op', 'sort_cond', '')
    field_default('op', 'sort_cond_vals', {})
    field_default('op', 'regions_to_plot', {})
    field_default('op', 'smooth_windowsize', 30)
    field_default('op', 'smooth_method', 'gaussian')
    field_default('op', 'trace_width', 1.5)
    field_default('op', 'epochs', [])
    field_default('op', 'epochs_to_label', [])
    field_default('op', 'epoch_colors', [])
    field_default('op', 'epoch_alpha', 0.12)
    field_default('op', 'epoch_label_height', 0.92)
    field_default('op', 'epoch_label_fontsize', 8)
    field_default('op', 'cmapname', 'jet')
    field_default('op', 'leg_pos_adjust', 0.22)
    field_default('op', 'y_ax_hardlims', [])
    field_default('op', 'samp_period', nan)
    
    %% ========== STEP 1: Filter electrodes ==========
    % Filter by responsiveness
    if op.analyze_responsive_elcs_only
         resp = resp(resp.rspv,:);
    end
    
    % Filter by tuning
    if op.analyze_tuned_elcs_only
        if ~isfield(op, 'tuning_param')
            error('op.analyze_tuned_elcs_only==1 but op.tuning_param not specified')
        end
        
        if ~ismember(op.tuning_param, resp.Properties.VariableNames)
            error('Tuning parameter "%s" not found in resp table', op.tuning_param)
        end
        
        % Handle nested indexing (e.g., {'col_name', index})
        if iscell(op.tuning_param)
            col_name = op.tuning_param{1};
            idx = op.tuning_param{2};
            tuned_elcs = resp{:, col_name}(:, idx) < op.tuning_alpha;
        else
            tuned_elcs = resp{:, op.tuning_param} < op.tuning_alpha;
        end
        
        resp = resp(tuned_elcs, :);
    end
    
    if height(resp) == 0
        error('No electrodes passed filtering criteria')
    end
    
    %% ========== STEP 2: Define brain regions ==========
    [resp, op] = define_brain_regions(resp, op);
    
    %% ========== STEP 3: Process each electrode ==========
    n_elc = height(resp);
    resp_by_cond = cell(n_elc, 1);
    times_aligned = cell(n_elc, 1);
    
    for ielc = 1:n_elc
        % Get this electrode's trials
        subind = find(string(subs.sub) == resp.sub{ielc});
        if isempty(subind)
            error('Subject %s not found in subs table', resp.sub{ielc})
        end
        trials_this_elc = subs.trials{subind};
        
        % Add warped timecourse as resp_aligned
        resp_tc = resp.timecourse{ielc};
        if iscell(resp_tc)
            trials_this_elc.resp_aligned = cell2mat(resp_tc);
        else
            trials_this_elc.resp_aligned = resp_tc;
        end
        
        % Extract times
        if ismember('times', resp.Properties.VariableNames)
            trials_this_elc.times = resp.times{ielc};
        end
        
        % Group by condition
        op_temp = op;
        op_temp.do_condition_sorting = 1;  % Always sort within electrode
        op_temp.sort_cond = op.sort_cond; 
        
        [~, align_stats, resp_grpd, ~] = sort_responses_by_condition_prealigned(trials_this_elc, op_temp);
        
        resp_by_cond{ielc} = resp_grpd;
        times_aligned{ielc} = align_stats.times_aligned;
        
    end
    
    %% ========== STEP 4: Aggregate by condition & electrode ==========
     cond_elc_data = [];
    for ielc = 1:n_elc
        resp_grpd_this = resp_by_cond{ielc};
        
        % Add metadata columns
        n_conds = height(resp_grpd_this);
        sub_col = repmat({resp.sub{ielc}}, n_conds, 1);
        chan_col = repmat(resp.chan(ielc), n_conds, 1);
        region_col = repmat({resp.region{ielc}}, n_conds, 1);
        times_col = repmat({times_aligned{ielc}}, n_conds, 1);
        
        metadata = table(sub_col, chan_col, region_col, times_col, ...
            'VariableNames', {'sub', 'chan', 'region', 'times'});
        
        cond_elc_data = [cond_elc_data; [metadata, resp_grpd_this]];
    end
    
    %% ========== STEP 5: Filter regions to plot if requested & Setup plotting ==========
    if ~isempty(op.regions_to_plot)
        if ~iscell(op.regions_to_plot)
            error('op.regions_to_plot must be a cell array of strings')
        end
        
        region_names_all = op.regiondef.region;
        valid_indices = [];
        for i = 1:length(op.regions_to_plot)
            req_reg = op.regions_to_plot{i};
            idx = find(strcmp(region_names_all, req_reg), 1);
            if isempty(op.regiondef.n_elcs(idx))
                fprintf(1, '    Warning: Specified region "%s" does not cover any electrodes in the resp table.\n', req_reg);
            else
                valid_indices(end+1) = idx;
            end
        end
        
        if isempty(valid_indices)
            error('None of the specified regions in op.regions_to_plot contain valid electrodes.')
        end
        
        op.regiondef = op.regiondef(valid_indices, :);
        op.nregions = height(op.regiondef);
    end

    if op.newfig
        hfig = figure('color', 'w', 'WindowState', 'maximized');
    end
    
    if op.analyze_tuned_elcs_only
        figtitle = sprintf('Electrodes tuned to %s (p<%g)', op.tuning_param, op.tuning_alpha);
    else
        figtitle = 'No electrode tuning criteria';
    end
    sgtitle(figtitle, 'FontSize', 14, 'FontWeight', 'bold')
    
    % Calculate subplot grid
    r = 1:op.nregions;
    c = ceil(op.nregions ./ r);
    [~, idx] = min(abs(c ./ r - 1/op.row_col_ratio));
    op.n_plot_rows = r(idx);
    op.n_plot_cols = c(idx);
    
    % Use tiledlayout for better spacing control
    tiledlayout(op.n_plot_rows, op.n_plot_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    %% ========== STEP 6: Plot by region ==========
   
    cond_elc_rgn = cell(op.nregions, 1);
    align_stats_rgn = cell(op.nregions, 1);
    resp_grpd_rgn = cell(op.nregions, 1);
    cfg_rgn = cell(op.nregions, 1);
    
    region_names = op.regiondef.region;
    
    for iregion = 1:op.nregions
        this_region = region_names{iregion};
        n_elcs_this_region = op.regiondef.n_elcs(iregion);
        
        % Create nexttile instead of subplot
        hax = nexttile;
        
        % Extract data for this region
        cond_elc_rgn{iregion} = cond_elc_data(strcmp(cond_elc_data.region, this_region), :);
        
        % Only plot if there's data AND enough electrodes
        if height(cond_elc_rgn{iregion}) > 0 && n_elcs_this_region >= 2
            
            % Aggregate responses by condition at regional level
            [align_stats_rgn{iregion}, resp_grpd_rgn{iregion}] = ...
                aggregate_regional_responses(cond_elc_rgn{iregion}, op);
            
            cfg_rgn{iregion} = op;
            cfg_rgn{iregion}.newfig = 0;
            cfg_rgn{iregion}.do_condition_sorting = 0;  % Already sorted
            cfg_rgn{iregion}.align_stats = align_stats_rgn{iregion};
            cfg_rgn{iregion}.resp_grpd = resp_grpd_rgn{iregion};
            
            % Get a trials table for plotting structure
            first_sub_idx = find(~cellfun(@isempty, subs.trials), 1);
            trials_for_plotting = subs.trials{first_sub_idx};
            
            plot_resp_timecourse(trials_for_plotting, cfg_rgn{iregion});
            
        else
            % Empty region - just add labels
            set(hax, 'XTickLabel', [])
            set(hax, 'YTickLabel', [])
            hold(hax, 'on')
            plot(hax, [0 1], [0 1], 'w')  % Invisible line to set up axes
            hold(hax, 'off')
        end
        
        % Add title for ALL regions (regardless of data)
        title_str = sprintf('%s (n=%d elc)', this_region, n_elcs_this_region);
        title(hax, title_str, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
        
    end
    
end

%% ========== HELPER: Sort by condition (pre-warped) ==========
function [trials_out, align_stats, resp_grpd, op_out] = ...
    sort_responses_by_condition_prealigned(trials, op)
    
    field_default('op', 'sort_cond_vals', {})
    field_default('op', 'samp_period', nan)
    field_default('op', 'sort_cond', '')
    
    %% Extract resp_aligned and convert to matrix if needed
    resp_aligned = trials.resp_aligned;
    
    % Convert cell array to matrix if needed
    if iscell(resp_aligned)
        resp_aligned = cell2mat(resp_aligned);
    end
    
    % Ensure it's a 2D matrix (trials x timepoints)
    if size(resp_aligned, 1) == 1
        resp_aligned = resp_aligned';
    end
    
    % Create minimal align_stats
    if ismember('times', trials.Properties.VariableNames)
        if iscell(trials.times)
            times_aligned = trials.times{1};
        else
            times_aligned = trials.times(1, :);
        end
    else
        % Estimate from response shape
        times_aligned = 1:size(resp_aligned, 2);
    end
    
    if isnan(op.samp_period)
        op.samp_period = mean(diff(times_aligned));
    end
    
    align_stats.times_aligned = times_aligned;
    align_stats.n_tpoints_pre_fixed = length(times_aligned);
    align_stats.n_tpoints_post_fixed = 0;
    align_stats.samp_period = op.samp_period;
    align_stats.mean = mean(resp_aligned, 1, 'omitnan');
    align_stats.std = std(resp_aligned, [], 'omitnan');
    align_stats.sem = align_stats.std ./ sqrt(sum(~isnan(resp_aligned), 1));
    align_stats.sem_lims = [align_stats.mean + align_stats.sem; 
                            align_stats.mean - align_stats.sem];
    
    %% Process sorting condition if specified
    if strcmp(op.sort_cond, '')
        % No sorting - return everything as single group
        resp_grpd = table();
        trials_out = trials;
        op_out = op;
        return
    end
    
    % Check that sort_cond exists
    if iscell(op.sort_cond)
        sort_col = op.sort_cond{1};
        sort_idx = op.sort_cond{2};
    else
        sort_col = op.sort_cond;
    end
    
    if ~ismember(sort_col, trials.Properties.VariableNames)
        error('Sort condition "%s" not found in trials table', sort_col)
    end
    
    % Extract sort condition values
    if iscell(op.sort_cond)
        trials.sort_cond = trials{:, sort_col}(:, sort_idx);
    else
        trials.sort_cond = trials{:, op.sort_cond};
    end
    
    % Get unique condition values
    if isempty(op.sort_cond_vals)
        op.sort_cond_vals = unique(trials.sort_cond);
    end
    
    % Convert numeric to string if needed
    if isnumeric(op.sort_cond_vals)
        op.sort_cond_vals(isnan(op.sort_cond_vals)) = [];
        op.sort_cond_vals = cellstr(string(op.sort_cond_vals));
    end
    
    [~, trial_cond_ind] = ismember(string(trials.sort_cond), op.sort_cond_vals);
    nconds = length(op.sort_cond_vals);
    
    %% Group responses by condition
    celcol = cell(nconds, 1);
    resp_grpd = table(reshape(op.sort_cond_vals, [], 1), celcol, celcol, ...
        'VariableNames', {'condval', 'resp', 'resp_mean'});
    
    for icond = 1:nconds
        these_trial_inds = trial_cond_ind == icond;
        resp_grpd.resp{icond} = resp_aligned(these_trial_inds, :);
        resp_grpd.resp_mean{icond} = mean(resp_grpd.resp{icond}, 1, 'omitnan');
        resp_grpd.std{icond} = std(resp_grpd.resp{icond}, 'omitnan');
        resp_grpd.n_good_trials{icond} = sum(~isnan(resp_grpd.resp{icond}));
        resp_grpd.sem{icond} = resp_grpd.std{icond} ./ sqrt(resp_grpd.n_good_trials{icond});
    end
    
    trials_out = trials;
    op_out = op;
    
end

%% ========== HELPER: Aggregate responses by condition at regional level ==========
function [align_stats, resp_grpd] = aggregate_regional_responses(cond_elc_data, op)
    % Aggregate electrode-level grouped responses into regional averages
    
    % Get unique conditions, respecting op.sort_cond_vals if provided
    if ~isempty(op.sort_cond_vals)
        if isnumeric(op.sort_cond_vals)
            requested_conds = cellstr(string(op.sort_cond_vals));
        else
            requested_conds = op.sort_cond_vals;
        end
        present_conds = unique(cond_elc_data.condval, 'stable');
        [~, loc] = ismember(cellstr(string(requested_conds)), cellstr(string(present_conds)));
        loc = loc(loc > 0);
        unique_conds = present_conds(loc);
    else
        unique_conds = unique(cond_elc_data.condval, 'stable');
    end
    
    nconds = length(unique_conds);
    
    % Initialize resp_grpd table
    celcol = cell(nconds, 1);
    resp_grpd = table(reshape(unique_conds, [], 1), celcol, celcol, ...
        'VariableNames', {'condval', 'resp', 'resp_mean'});
    
    % Get times from first electrode
    times_aligned = cond_elc_data.times{1};
    
    % Convert times to trial-relative (start at 0)
    times_aligned = times_aligned - times_aligned(1);
    
    % Aggregate each condition
    for icond = 1:nconds
        cond_name = unique_conds{icond};
        cond_mask = strcmp(cond_elc_data.condval, cond_name);
        cond_data = cond_elc_data(cond_mask, :);
        
        % Stack all electrode means for this condition
        resp_means = cell2mat(cond_data.resp_mean);
        
        % Compute regional statistics using nanmean
        resp_grpd.resp{icond} = resp_means;
        resp_grpd.resp_mean{icond} = nanmean(resp_means, 1);
        resp_grpd.std{icond} = nanstd(resp_means);
        resp_grpd.n_good_trials{icond} = sum(~isnan(resp_means));
        resp_grpd.sem{icond} = resp_grpd.std{icond} ./ sqrt(resp_grpd.n_good_trials{icond});
    end
    
    % Create align_stats
    if isnan(op.samp_period)
        op.samp_period = mean(diff(times_aligned));
    end
    
    align_stats.times_aligned = times_aligned;
    align_stats.n_tpoints_pre_fixed = length(times_aligned);
    align_stats.n_tpoints_post_fixed = 0;
    align_stats.samp_period = op.samp_period;
    
    % Use nanmean for grand average
    all_means = cell2mat(resp_grpd.resp_mean);
    align_stats.mean = nanmean(all_means, 1);
    align_stats.std = nanstd(all_means);
    align_stats.sem = align_stats.std ./ sqrt(nconds);
    align_stats.sem_lims = [align_stats.mean + align_stats.sem;
                            align_stats.mean - align_stats.sem];
    
end