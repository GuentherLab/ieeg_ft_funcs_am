function [cond_elc_rgn, align_stats_rgn, resp_grpd_rgn, cfg_rgn] = ...
    combine_plot_electrode_timecourses(resp, subs, op)
    
    %% COMBINE_PLOT_ELECTRODE_TIMECOURSES
    % Main function for visualizing pre-warped electrode timecourses
    
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
    field_default('op', 'y_ax_hardlims', [])
    field_default('op', 'samp_period', nan)
    field_default('op', 'min_elcs_per_region', 1)
    
    %% Extract Master Reference Time Vector
    % Find a valid seconds-based time vector across electrodes to use as a template
    master_times = [];
    if ismember('times', resp.Properties.VariableNames)
        for i = 1:height(resp)
            t_cand = resp.times{i};
            if ~isempty(t_cand) && isnumeric(t_cand) && max(t_cand) < 50
                master_times = t_cand;
                break;
            end
        end
    end
    
    if isnan(op.samp_period) && ~isempty(master_times) && length(master_times) > 1
        op.samp_period = mean(diff(master_times));
    end
    op.master_times = master_times;

    %% ========== STEP 1: Filter electrodes ==========
    if op.analyze_responsive_elcs_only
         resp = resp(resp.rspv, :);
    end
    
    if op.analyze_tuned_elcs_only
        if ~isfield(op, 'tuning_param')
            error('op.analyze_tuned_elcs_only==1 but op.tuning_param not specified')
        end
        
        if ~ismember(op.tuning_param, resp.Properties.VariableNames)
            error('Tuning parameter "%s" not found in resp table', op.tuning_param)
        end
        
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
        subind = find(string(subs.sub) == string(resp.sub{ielc}));
        if isempty(subind)
            error('Subject %s not found in subs table', resp.sub{ielc})
        end
        trials_this_elc = subs.trials{subind};
        
        resp_tc = resp.timecourse{ielc};
        if iscell(resp_tc)
            if isrow(resp_tc)
                resp_tc = resp_tc';
            end
            trials_this_elc.resp_aligned = cell2mat(resp_tc);
        else
            trials_this_elc.resp_aligned = resp_tc;
        end
        
        if ismember('times', resp.Properties.VariableNames) && ~isempty(resp.times{ielc})
            trials_this_elc.times = resp.times{ielc};
        end
        
        op_temp = op;
        op_temp.do_condition_sorting = 1;
        op_temp.sort_cond = op.sort_cond; 
        
        [~, align_stats, resp_grpd, ~] = sort_responses_by_condition_prealigned(trials_this_elc, op_temp);
        
        resp_by_cond{ielc} = resp_grpd;
        times_aligned{ielc} = align_stats.times_aligned;
    end
    
    %% ========== STEP 4: Aggregate by condition & electrode ==========
    cond_elc_data = [];
    for ielc = 1:n_elc
        resp_grpd_this = resp_by_cond{ielc};
        
        n_conds = height(resp_grpd_this);
        sub_col = repmat({resp.sub{ielc}}, n_conds, 1);
        chan_col = repmat(resp.chan(ielc), n_conds, 1);
        region_col = repmat({resp.region{ielc}}, n_conds, 1);
        times_col = repmat({times_aligned{ielc}}, n_conds, 1);
        
        metadata = table(sub_col, chan_col, region_col, times_col, ...
            'VariableNames', {'sub', 'chan', 'region', 'times'});
        
        cond_elc_data = [cond_elc_data; [metadata, resp_grpd_this]];
    end
    
    %% ========== STEP 5: Filter regions & Setup layout ==========
    if ~isempty(op.regions_to_plot)
        if ~iscell(op.regions_to_plot)
            error('op.regions_to_plot must be a cell array of strings')
        end
        
        region_names_all = op.regiondef.region;
        valid_indices = [];
        for i = 1:length(op.regions_to_plot)
            req_reg = op.regions_to_plot{i};
            idx = find(strcmp(region_names_all, req_reg), 1);
            if isempty(idx) || isempty(op.regiondef.n_elcs(idx)) || op.regiondef.n_elcs(idx) == 0
                fprintf(1, '    Warning: Specified region "%s" does not cover any electrodes.\n', req_reg);
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
    
    r = 1:op.nregions;
    c = ceil(op.nregions ./ r);
    [~, idx] = min(abs(c ./ r - 1/op.row_col_ratio));
    op.n_plot_rows = r(idx);
    op.n_plot_cols = c(idx);
    
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
        
        hax = nexttile;
        cond_elc_rgn{iregion} = cond_elc_data(strcmp(cond_elc_data.region, this_region), :);
        
        if height(cond_elc_rgn{iregion}) > 0 && n_elcs_this_region >= op.min_elcs_per_region
            
            [align_stats_rgn{iregion}, resp_grpd_rgn{iregion}] = ...
                aggregate_regional_responses(cond_elc_rgn{iregion}, op);
            
            cfg_rgn{iregion} = op;
            cfg_rgn{iregion}.newfig = 0;
            cfg_rgn{iregion}.do_condition_sorting = 0;
            cfg_rgn{iregion}.align_stats = align_stats_rgn{iregion};
            cfg_rgn{iregion}.resp_grpd = resp_grpd_rgn{iregion};
            
            first_sub_idx = find(~cellfun(@isempty, subs.trials), 1);
            trials_for_plotting = subs.trials{first_sub_idx};
            
            plot_resp_timecourse(trials_for_plotting, cfg_rgn{iregion});
            
        else
            set(hax, 'XTickLabel', [])
            set(hax, 'YTickLabel', [])
            hold(hax, 'on')
            plot(hax, [0 1], [0 1], 'w')
            hold(hax, 'off')
        end
        
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
    
    resp_aligned = trials.resp_aligned;
    if iscell(resp_aligned)
        if isrow(resp_aligned)
            resp_aligned = resp_aligned';
        end
        resp_aligned = cell2mat(resp_aligned);
    end
    
    n_pts = size(resp_aligned, 2);
    
    %% Extract and validate time vector
    times_aligned = [];
    if ismember('times', trials.Properties.VariableNames) && ~isempty(trials.times)
        if iscell(trials.times) && ~isempty(trials.times{1})
            times_aligned = trials.times{1};
        elseif isnumeric(trials.times)
            times_aligned = trials.times(1, :);
        end
    end
    
    % Fallback 1: Master template
    if (isempty(times_aligned) || max(times_aligned) > 50) && isfield(op, 'master_times') && ~isempty(op.master_times) && length(op.master_times) == n_pts
        times_aligned = op.master_times;
    end
    
    % Fallback 2: Calculate from samp_period
    if isempty(times_aligned)
        if ~isnan(op.samp_period) && op.samp_period < 0.5
            times_aligned = (0:n_pts-1) * op.samp_period;
        else
            times_aligned = 1:n_pts;
        end
    end
    
    % Ensure times_aligned is in seconds (not sample indices 1..650)
    if max(times_aligned) > 50 && ~isnan(op.samp_period) && op.samp_period < 0.5
        times_aligned = (0:n_pts-1) * op.samp_period;
    end
    
    times_aligned = times_aligned - times_aligned(1);
    
    if isnan(op.samp_period)
        op.samp_period = mean(diff(times_aligned));
    end
    
    align_stats.times_aligned = times_aligned;
    align_stats.n_tpoints_pre_fixed = length(times_aligned);
    align_stats.n_tpoints_post_fixed = 0;
    align_stats.samp_period = op.samp_period;
    align_stats.mean = nanmean(resp_aligned, 1);
    align_stats.std = nanstd(resp_aligned);
    align_stats.sem = align_stats.std ./ sqrt(sum(~isnan(resp_aligned), 1));
    align_stats.sem_lims = [align_stats.mean + align_stats.sem; 
                            align_stats.mean - align_stats.sem];
    
    if isempty(op.sort_cond) || strcmp(op.sort_cond, '')
        celcol = cell(1, 1);
        resp_grpd = table({'All'}, celcol, celcol, ...
            'VariableNames', {'condval', 'resp', 'resp_mean'});
        resp_grpd.resp{1} = resp_aligned;
        resp_grpd.resp_mean{1} = nanmean(resp_aligned, 1);
        resp_grpd.std{1} = nanstd(resp_aligned);
        resp_grpd.n_good_trials{1} = sum(~isnan(resp_aligned), 1);
        resp_grpd.sem{1} = resp_grpd.std{1} ./ sqrt(resp_grpd.n_good_trials{1});
        
        trials_out = trials;
        op_out = op;
        return
    end
    
    if iscell(op.sort_cond)
        sort_col = op.sort_cond{1};
        sort_idx = op.sort_cond{2};
    else
        sort_col = op.sort_cond;
    end
    
    if ~ismember(sort_col, trials.Properties.VariableNames)
        error('Sort condition "%s" not found in trials table', sort_col)
    end
    
    if iscell(op.sort_cond)
        trials.sort_cond = trials{:, sort_col}(:, sort_idx);
    else
        trials.sort_cond = trials{:, sort_col};
    end
    
    if ~isempty(op.sort_cond_vals)
        cond_vals_str = cellstr(string(op.sort_cond_vals));
        [keep_trials, trial_cond_ind] = ismember(string(trials.sort_cond), cond_vals_str);
        
        trials = trials(keep_trials, :);
        resp_aligned = resp_aligned(keep_trials, :);
        trial_cond_ind = trial_cond_ind(keep_trials);
        
        op.sort_cond_vals = cond_vals_str;
    else
        unique_vals = unique(trials.sort_cond);
        if isnumeric(unique_vals)
            unique_vals = unique_vals(~isnan(unique_vals));
        end
        op.sort_cond_vals = cellstr(string(unique_vals));
        
        [~, trial_cond_ind] = ismember(string(trials.sort_cond), string(op.sort_cond_vals));
    end
    
    nconds = length(op.sort_cond_vals);
    
    celcol = cell(nconds, 1);
    resp_grpd = table(reshape(op.sort_cond_vals, [], 1), celcol, celcol, ...
        'VariableNames', {'condval', 'resp', 'resp_mean'});
    
    for icond = 1:nconds
        these_trial_inds = trial_cond_ind == icond;
        
        if nnz(these_trial_inds) > 0
            resp_grpd.resp{icond} = resp_aligned(these_trial_inds, :);
            resp_grpd.resp_mean{icond} = nanmean(resp_grpd.resp{icond}, 1);
            resp_grpd.std{icond} = nanstd(resp_grpd.resp{icond});
            resp_grpd.n_good_trials{icond} = sum(~isnan(resp_grpd.resp{icond}));
            resp_grpd.sem{icond} = resp_grpd.std{icond} ./ sqrt(resp_grpd.n_good_trials{icond});
        else
            resp_grpd.resp{icond} = [];
            resp_grpd.resp_mean{icond} = [];
            resp_grpd.std{icond} = [];
            resp_grpd.n_good_trials{icond} = 0;
            resp_grpd.sem{icond} = [];
        end
    end
    
    trials_out = trials;
    op_out = op;
end

%% ========== HELPER: Aggregate responses by condition at regional level ==========
function [align_stats, resp_grpd] = aggregate_regional_responses(cond_elc_data, op)
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
    
    celcol = cell(nconds, 1);
    resp_grpd = table(reshape(unique_conds, [], 1), celcol, celcol, ...
        'VariableNames', {'condval', 'resp', 'resp_mean'});
    
    % Find valid time vector across electrodes in this region
    times_aligned = [];
    for irow = 1:height(cond_elc_data)
        t_cand = cond_elc_data.times{irow};
        if ~isempty(t_cand) && max(t_cand) < 50
            times_aligned = t_cand;
            break;
        end
    end
    
    % Use master_times fallback if necessary
    if isempty(times_aligned) && isfield(op, 'master_times') && ~isempty(op.master_times)
        times_aligned = op.master_times;
    end
    
    % Use length of first response vector as last fallback
    if isempty(times_aligned)
        for irow = 1:height(cond_elc_data)
            if ~isempty(cond_elc_data.resp_mean{irow})
                n_pts = length(cond_elc_data.resp_mean{irow});
                if ~isnan(op.samp_period) && op.samp_period < 0.5
                    times_aligned = (0:n_pts-1) * op.samp_period;
                else
                    times_aligned = 1:n_pts;
                end
                break;
            end
        end
    end
    
    % Convert to relative time starting at 0
    times_aligned = times_aligned - times_aligned(1);
    
    % Rescale sample indices to seconds if needed
    if max(times_aligned) > 50 && ~isnan(op.samp_period) && op.samp_period < 0.5
        times_aligned = (0:length(times_aligned)-1) * op.samp_period;
    end
    
    for icond = 1:nconds
        cond_name = unique_conds{icond};
        cond_mask = strcmp(cond_elc_data.condval, cond_name);
        cond_data = cond_elc_data(cond_mask, :);
        
        resp_means = cell2mat(cond_data.resp_mean);
        
        resp_grpd.resp{icond} = resp_means;
        resp_grpd.resp_mean{icond} = nanmean(resp_means, 1);
        resp_grpd.std{icond} = nanstd(resp_means);
        resp_grpd.n_good_trials{icond} = sum(~isnan(resp_means));
        resp_grpd.sem{icond} = resp_grpd.std{icond} ./ sqrt(resp_grpd.n_good_trials{icond});
    end
    
    if isnan(op.samp_period)
        op.samp_period = mean(diff(times_aligned));
    end
    
    align_stats.times_aligned = times_aligned;
    align_stats.n_tpoints_pre_fixed = length(times_aligned);
    align_stats.n_tpoints_post_fixed = 0;
    align_stats.samp_period = op.samp_period;
    
    all_means = cell2mat(resp_grpd.resp_mean);
    align_stats.mean = nanmean(all_means, 1);
    align_stats.std = nanstd(all_means);
    align_stats.sem = align_stats.std ./ sqrt(nconds);
    align_stats.sem_lims = [align_stats.mean + align_stats.sem;
                            align_stats.mean - align_stats.sem];
end