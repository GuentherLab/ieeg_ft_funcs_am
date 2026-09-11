function [electrode_plot_data, align_stats_elc, resp_grpd_elc, cfg_elc, fig_handles, ax_handles] = ...
    plot_electrode_timecourses(resp, subs, op)
% PLOT_ELECTRODE_TIMECOURSES
%
% Plot pre-warped electrode response timecourses with ONE ELECTRODE PER
% SUBPLOT.

    %% --------------------------------------------------------------------
    % Defaults
    % ---------------------------------------------------------------------

    if nargin < 3 || isempty(op)
        op = struct();
    end

    op = set_default(op, 'newfig', true);

    op = set_default(op, 'analyze_responsive_elcs_only', true);
    op = set_default(op, 'analyze_tuned_elcs_only', true);
    op = set_default(op, 'tuning_alpha', 0.05);

    op = set_default(op, 'subjects_to_plot', {});
    op = set_default(op, 'regions_to_plot', {});
    op = set_default(op, 'min_elcs_per_region', 1);

    op = set_default(op, 'sort_cond', '');
    op = set_default(op, 'sort_cond_vals', {});
    op = set_default(op, 'condval_inds_to_plot', []);

    op = set_default(op, 'plotrows', 4);
    op = set_default(op, 'plotcolumns', 4);

    op = set_default(op, 'screen_order', []);

    op = set_default(op, 'smooth_timecourses', true);
    op = set_default(op, 'smooth_windowsize', 30);
    op = set_default(op, 'smooth_method', 'gaussian');

    op = set_default(op, 'trace_width', 1.5);
    op = set_default(op, 'cmapname', 'jet');

    op = set_default(op, 'epochs', []);
    op = set_default(op, 'epochs_to_label', []);
    op = set_default(op, 'epoch_colors', []);
    op = set_default(op, 'epoch_alpha', 0.12);
    op = set_default(op, 'epoch_label_height', 0.92);
    op = set_default(op, 'epoch_label_fontsize', 8);
    op = set_default(op, 'epoch_divider_linestyle', ':');
    op = set_default(op, 'epoch_divider_linewidth', 1.5);

    op = set_default(op, 'y_ax_hardlims', []);
    op = set_default(op, 'yline_zero_width', 0.25);
    op = set_default(op, 'yline_zero_color', [0.8 0.8 0.8]);
    op = set_default(op, 'yline_zero_style', '-');

    op = set_default(op, 'samp_period', nan);

    op = set_default(op, 'sem_fill_color', [0.35 0.35 0.35]);
    op = set_default(op, 'sem_fill_alpha', 0.34);

    op = set_default(op, 'show_legend', true);
    op = set_default(op, 'legend_location', 'northeast');
    op = set_default(op, 'legend_fontsize', 6);
    op = set_default(op, 'legend_box', 'on');
    op = set_default(op, 'legend_auto_avoid_overlap', true);

    if op.plotrows < 1 || op.plotcolumns < 1
        error('op.plotrows and op.plotcolumns must both be positive integers.')
    end

    op.plotrows = round(op.plotrows);
    op.plotcolumns = round(op.plotcolumns);

    %% --------------------------------------------------------------------
    % Extract a seconds-based master time vector, if one exists
    % ---------------------------------------------------------------------

    master_times = [];

    if ismember('times', resp.Properties.VariableNames)
        for irow = 1:height(resp)
            time_candidate = extract_resp_time_vector(resp, irow);

            if is_plausible_seconds_time_vector(time_candidate)
                master_times = time_candidate(:)';
                break
            end
        end
    end

    if isnan(op.samp_period) && ~isempty(master_times) && numel(master_times) > 1
        op.samp_period = mean(diff(master_times));
    end

    op.master_times = master_times;

    %% --------------------------------------------------------------------
    % STEP 1: Electrode-level filtering
    % ---------------------------------------------------------------------

    electrode_keep_mask = true(height(resp), 1);

    if op.analyze_responsive_elcs_only
        if ~ismember('rspv', resp.Properties.VariableNames)
            error('op.analyze_responsive_elcs_only is true, but resp.rspv was not found.')
        end

        electrode_keep_mask = electrode_keep_mask & logical(resp.rspv);
    end

    if op.analyze_tuned_elcs_only
        if ~isfield(op, 'tuning_param') || isempty(op.tuning_param)
            error('op.analyze_tuned_elcs_only is true, but op.tuning_param was not specified.')
        end

        if iscell(op.tuning_param)
            tuning_column_name = op.tuning_param{1};
            tuning_column_index = op.tuning_param{2};
        else
            tuning_column_name = char(string(op.tuning_param));
            tuning_column_index = [];
        end

        if ~ismember(tuning_column_name, resp.Properties.VariableNames)
            error('Tuning parameter "%s" was not found in resp table.', tuning_column_name)
        end

        tuning_values = resp{:, tuning_column_name};

        if iscell(tuning_values)
            tuning_values = cell2mat(tuning_values);
        end

        if isempty(tuning_column_index)
            if size(tuning_values, 2) ~= 1
                error(['Tuning parameter "%s" has multiple columns. ', ...
                       'Use op.tuning_param = {''%s'', column_index}.'], ...
                       tuning_column_name, tuning_column_name)
            end

            tuned_electrode_mask = tuning_values < op.tuning_alpha;
        else
            tuned_electrode_mask = tuning_values(:, tuning_column_index) < op.tuning_alpha;
        end

        electrode_keep_mask = electrode_keep_mask & tuned_electrode_mask;
    end

    if ~isempty(op.subjects_to_plot)
        if ~ismember('sub', resp.Properties.VariableNames)
            error('op.subjects_to_plot was specified, but resp.sub was not found.')
        end

        requested_subjects = safe_string_column(op.subjects_to_plot);
        response_subjects = safe_string_column(resp.sub);

        electrode_keep_mask = electrode_keep_mask & ismember(response_subjects, requested_subjects);
    end

    resp = resp(electrode_keep_mask, :);

    if height(resp) == 0
        error('No electrodes passed responsive/tuning/subject filtering criteria.')
    end

    %% --------------------------------------------------------------------
    % STEP 2: Define or verify brain regions
    % ---------------------------------------------------------------------

    if exist('define_brain_regions', 'file') == 2
        [resp, op] = define_brain_regions(resp, op);
    else
        if ~ismember('region', resp.Properties.VariableNames)
            error(['define_brain_regions.m was not found, and resp.region does not exist. ', ...
                   'Either add define_brain_regions.m to the path or provide resp.region.'])
        end

        op.regiondef = make_region_definition_table(resp, {});
        op.nregions = height(op.regiondef);
    end

    %% --------------------------------------------------------------------
    % STEP 3: Region inclusion filter
    % ---------------------------------------------------------------------

    if ~isempty(op.regions_to_plot)
        if ~iscell(op.regions_to_plot) && ~isstring(op.regions_to_plot)
            error('op.regions_to_plot must be a cell array of strings or a string array.')
        end

        requested_regions = safe_string_column(op.regions_to_plot);
        available_regions = safe_string_column(resp.region);

        missing_regions = setdiff(requested_regions, unique(available_regions));

        for imissing = 1:numel(missing_regions)
            fprintf(1, 'Warning: Specified region "%s" does not cover any electrodes after filtering.\n', ...
                missing_regions(imissing));
        end

        region_keep_mask = ismember(available_regions, requested_regions);
        resp = resp(region_keep_mask, :);

        if height(resp) == 0
            error('None of the specified regions in op.regions_to_plot contain valid electrodes.')
        end
    end

    if op.min_elcs_per_region > 1
        region_names_after_filtering = safe_string_column(resp.region);
        unique_region_names = unique(region_names_after_filtering, 'stable');
        unique_region_names = unique_region_names(unique_region_names ~= "");

        region_has_enough_electrodes = false(size(region_names_after_filtering));

        for iregion = 1:numel(unique_region_names)
            this_region_name = unique_region_names(iregion);
            this_region_mask = region_names_after_filtering == this_region_name;
            this_region_count = nnz(this_region_mask);

            if this_region_count >= op.min_elcs_per_region
                region_has_enough_electrodes = region_has_enough_electrodes | this_region_mask;
            end
        end

        resp = resp(region_has_enough_electrodes, :);

        if height(resp) == 0
            error('No electrodes remain after applying op.min_elcs_per_region.')
        end
    end

    op.regiondef = make_region_definition_table(resp, op.regions_to_plot);
    op.nregions = height(op.regiondef);

    %% --------------------------------------------------------------------
    % STEP 4: Build trial-level grouped data for each electrode
    % ---------------------------------------------------------------------

    n_electrodes_to_plot = height(resp);

    trials_elc = cell(n_electrodes_to_plot, 1);
    align_stats_elc = cell(n_electrodes_to_plot, 1);
    resp_grpd_elc = cell(n_electrodes_to_plot, 1);
    cfg_elc = cell(n_electrodes_to_plot, 1);

    response_subjects = safe_string_column(resp.sub);
    subs_subjects = safe_string_column(subs.sub);

    for ielc = 1:n_electrodes_to_plot

        subject_id = response_subjects(ielc);
        subject_row_index = find(subs_subjects == subject_id, 1);

        if isempty(subject_row_index)
            error('Subject "%s" from resp.sub was not found in subs.sub.', subject_id)
        end

        trials_this_electrode = subs.trials{subject_row_index};

        if isempty(trials_this_electrode) || height(trials_this_electrode) == 0
            error('Subject "%s" has an empty trials table.', subject_id)
        end

        response_matrix = normalize_electrode_timecourse(resp.timecourse{ielc}, height(trials_this_electrode));

        if size(response_matrix, 1) ~= height(trials_this_electrode)
            error(['Electrode %d has %d response rows, but subject "%s" has %d trials. ', ...
                   'The electrode timecourse must have one row per trial.'], ...
                   ielc, size(response_matrix, 1), subject_id, height(trials_this_electrode))
        end

        trials_this_electrode.resp_aligned = response_matrix;

        electrode_time_vector = extract_resp_time_vector(resp, ielc);

        if ~isempty(electrode_time_vector)
            trials_this_electrode.times = repmat({electrode_time_vector(:)'}, height(trials_this_electrode), 1);
        end

        cfg_this_electrode = op;
        cfg_this_electrode.do_condition_sorting = true;

        [trials_this_electrode, align_stats_this_electrode, resp_grpd_this_electrode, cfg_this_electrode] = ...
            sort_responses_by_condition_prealigned_single_electrode(trials_this_electrode, cfg_this_electrode);

        trials_elc{ielc} = trials_this_electrode;
        align_stats_elc{ielc} = align_stats_this_electrode;
        resp_grpd_elc{ielc} = resp_grpd_this_electrode;

        cfg_this_electrode.newfig = false;
        cfg_this_electrode.do_condition_sorting = false;
        cfg_this_electrode.align_stats = align_stats_this_electrode;
        cfg_this_electrode.resp_grpd = resp_grpd_this_electrode;

        cfg_elc{ielc} = cfg_this_electrode;
    end

    %% --------------------------------------------------------------------
    % STEP 5: Create figures and plot each electrode
    % ---------------------------------------------------------------------

    n_tiles_per_figure = op.plotrows * op.plotcolumns;
    n_figures_needed = ceil(n_electrodes_to_plot / n_tiles_per_figure);

    fig_handles = gobjects(n_figures_needed, 1);
    ax_handles = gobjects(n_electrodes_to_plot, 1);

    for ifig = 1:n_figures_needed

        electrode_start_index = (ifig - 1) * n_tiles_per_figure + 1;
        electrode_stop_index = min(ifig * n_tiles_per_figure, n_electrodes_to_plot);

        fig_handles(ifig) = create_positioned_figure(ifig, op);

        figure(fig_handles(ifig));

        tile_layout_handle = tiledlayout(op.plotrows, op.plotcolumns, ...
            'Padding', 'compact', ...
            'TileSpacing', 'compact');

        figure_title = sprintf('Electrode timecourses: electrodes %d-%d of %d', ...
            electrode_start_index, electrode_stop_index, n_electrodes_to_plot);

        if op.analyze_tuned_elcs_only
            figure_title = sprintf('%s | tuned to %s, p < %g', ...
                figure_title, option_value_to_title_string(op.tuning_param), op.tuning_alpha);
        else
            figure_title = sprintf('%s | no electrode tuning criteria', figure_title);
        end

        sgtitle(tile_layout_handle, figure_title, ...
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');

        for ielc = electrode_start_index:electrode_stop_index

            tile_index_within_figure = ielc - electrode_start_index + 1;

            ax_handles(ielc) = nexttile(tile_layout_handle, tile_index_within_figure);

            axes(ax_handles(ielc)); %#ok<LAXES>

            plot_single_electrode_timecourse(trials_elc{ielc}, cfg_elc{ielc});

            electrode_title = make_electrode_title(resp, ielc);

            title(ax_handles(ielc), electrode_title, ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'Interpreter', 'none');
        end
    end

    %% --------------------------------------------------------------------
    % STEP 6: Return filtered resp table with plot bookkeeping
    % ---------------------------------------------------------------------

    electrode_plot_data = resp;
    electrode_plot_data.plot_order = (1:n_electrodes_to_plot)';
    electrode_plot_data.figure_index = ceil(electrode_plot_data.plot_order ./ n_tiles_per_figure);
    electrode_plot_data.tile_index = mod(electrode_plot_data.plot_order - 1, n_tiles_per_figure) + 1;
end

%% ========================================================================
% Helper: sort/group responses for one pre-warped electrode
% ========================================================================

function [trials_out, align_stats, resp_grpd, op_out] = ...
    sort_responses_by_condition_prealigned_single_electrode(trials, op)

    op = set_default(op, 'sort_cond_vals', {});
    op = set_default(op, 'samp_period', nan);
    op = set_default(op, 'sort_cond', '');

    resp_aligned = trials.resp_aligned;

    if iscell(resp_aligned)
        resp_aligned = cell2mat(resp_aligned);
    end

    n_timepoints = size(resp_aligned, 2);

    raw_times_aligned = extract_trial_time_vector(trials);

    [times_aligned, op] = standardize_time_vector_to_seconds(raw_times_aligned, n_timepoints, op);

    align_stats.times_aligned = times_aligned;
    align_stats.n_tpoints_pre_fixed = numel(times_aligned);
    align_stats.n_tpoints_post_fixed = 0;
    align_stats.samp_period = op.samp_period;

    align_stats.mean = nanmean(resp_aligned, 1);
    align_stats.std = nanstd(resp_aligned, [], 1);
    align_stats.n_good_trials = sum(~isnan(resp_aligned), 1);
    align_stats.sem = align_stats.std ./ sqrt(align_stats.n_good_trials);

    align_stats.sem_lims = [align_stats.mean + align_stats.sem;
                            align_stats.mean - align_stats.sem];

    if is_empty_sort_condition(op.sort_cond)
        resp_grpd = make_single_condition_response_group('All Trials', resp_aligned, align_stats);

        trials_out = trials;
        op_out = op;
        return
    end

    if iscell(op.sort_cond)
        sort_column_name = op.sort_cond{1};
        sort_column_index = op.sort_cond{2};
    else
        sort_column_name = char(string(op.sort_cond));
        sort_column_index = [];
    end

    if ~ismember(sort_column_name, trials.Properties.VariableNames)
        error('Sort condition "%s" was not found in trials table.', sort_column_name)
    end

    if isempty(sort_column_index)
        trials.sort_cond = trials{:, sort_column_name};
    else
        trials.sort_cond = trials{:, sort_column_name}(:, sort_column_index);
    end

    if ~isempty(op.sort_cond_vals)
        condition_values_to_plot = cellstr(safe_string_column(op.sort_cond_vals));

        [keep_trials, trial_condition_index] = ismember( ...
            safe_string_column(trials.sort_cond), ...
            string(condition_values_to_plot));

        trials = trials(keep_trials, :);
        resp_aligned = resp_aligned(keep_trials, :);
        trial_condition_index = trial_condition_index(keep_trials);

        op.sort_cond_vals = condition_values_to_plot;
    else
        unique_condition_values = unique(trials.sort_cond, 'stable');

        if isnumeric(unique_condition_values)
            unique_condition_values = unique_condition_values(~isnan(unique_condition_values));
        end

        op.sort_cond_vals = cellstr(safe_string_column(unique_condition_values));

        [~, trial_condition_index] = ismember( ...
            safe_string_column(trials.sort_cond), ...
            string(op.sort_cond_vals));
    end

    n_conditions = numel(op.sort_cond_vals);

    empty_cell_column = cell(n_conditions, 1);

    resp_grpd = table(reshape(op.sort_cond_vals, [], 1), ...
                      empty_cell_column, ...
                      empty_cell_column, ...
                      'VariableNames', {'condval', 'resp', 'resp_mean'});

    for icond = 1:n_conditions
        these_trial_indices = trial_condition_index == icond;

        if nnz(these_trial_indices) > 0
            responses_this_condition = resp_aligned(these_trial_indices, :);

            resp_grpd.resp{icond} = responses_this_condition;
            resp_grpd.resp_mean{icond} = nanmean(responses_this_condition, 1);
            resp_grpd.std{icond} = nanstd(responses_this_condition, [], 1);
            resp_grpd.n_good_trials{icond} = sum(~isnan(responses_this_condition), 1);
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

%% ========================================================================
% Helper: robustly standardize time vector to seconds
% ========================================================================

function [times_seconds, op] = standardize_time_vector_to_seconds(raw_times, n_timepoints, op)

    if nargin < 3 || isempty(op)
        op = struct();
    end

    op = set_default(op, 'samp_period', nan);
    op = set_default(op, 'master_times', []);
    op = set_default(op, 'epochs', []);

    raw_times = raw_times(:)';

    raw_is_usable_length = ~isempty(raw_times) && numel(raw_times) == n_timepoints;
    raw_is_seconds = raw_is_usable_length && is_plausible_seconds_time_vector(raw_times);
    raw_looks_like_samples = raw_is_usable_length && looks_like_sample_indices(raw_times);

    if raw_is_seconds
        times_seconds = raw_times - raw_times(1);

        if isnan(op.samp_period) && numel(times_seconds) > 1
            op.samp_period = mean(diff(times_seconds));
        end

        return
    end

    if isfield(op, 'master_times') && ...
            ~isempty(op.master_times) && ...
            numel(op.master_times) == n_timepoints && ...
            is_plausible_seconds_time_vector(op.master_times)

        times_seconds = op.master_times(:)' - op.master_times(1);

        if isnan(op.samp_period) && numel(times_seconds) > 1
            op.samp_period = mean(diff(times_seconds));
        end

        return
    end

    if ~isnan(op.samp_period) && op.samp_period > 0 && op.samp_period < 0.5
        times_seconds = (0:n_timepoints-1) * op.samp_period;
        return
    end

    inferred_sample_period = infer_sample_period_from_epochs(op, n_timepoints);

    if ~isnan(inferred_sample_period)
        op.samp_period = inferred_sample_period;
        times_seconds = (0:n_timepoints-1) * op.samp_period;
        return
    end

    if raw_is_usable_length
        times_seconds = raw_times - raw_times(1);

        if raw_looks_like_samples || max(times_seconds) > 50
            warning(['Time vector appears to be sample indices, but no valid ', ...
                     'op.samp_period, master_times, or op.epochs duration was available. ', ...
                     'Plot x-axis may be in samples, not seconds.'])
        end

        return
    end

    times_seconds = 1:n_timepoints;
end

function tf = is_plausible_seconds_time_vector(time_vector)

    tf = false;

    if isempty(time_vector) || ~isnumeric(time_vector)
        return
    end

    time_vector = time_vector(:)';

    if numel(time_vector) < 2
        return
    end

    if any(~isfinite(time_vector))
        return
    end

    time_range = max(time_vector) - min(time_vector);

    if time_range <= 0
        return
    end

    if max(abs(time_vector)) < 50 && time_range < 50
        tf = true;
    end
end

function tf = looks_like_sample_indices(time_vector)

    tf = false;

    if isempty(time_vector) || ~isnumeric(time_vector) || numel(time_vector) < 2
        return
    end

    time_vector = time_vector(:)';

    if any(~isfinite(time_vector))
        return
    end

    diffs = diff(time_vector);

    if isempty(diffs)
        return
    end

    median_step = median(diffs);

    if abs(median_step - 1) < 1e-6 && max(time_vector) > 50
        tf = true;
    end
end

function inferred_sample_period = infer_sample_period_from_epochs(op, n_timepoints)

    inferred_sample_period = nan;

    if ~isfield(op, 'epochs') || isempty(op.epochs)
        return
    end

    if ~istable(op.epochs) || ~ismember('dur_fix', op.epochs.Properties.VariableNames)
        return
    end

    epoch_durations = op.epochs.dur_fix(:);
    total_epoch_duration = sum(epoch_durations);

    if isempty(total_epoch_duration) || ~isfinite(total_epoch_duration) || total_epoch_duration <= 0
        return
    end

    if n_timepoints > 1
        inferred_sample_period = total_epoch_duration / (n_timepoints - 1);
    end
end

%% ========================================================================
% Helper: plot one electrode into current axes
% ========================================================================

function [trials_out, resp_grpd, align_stats, op_out] = plot_single_electrode_timecourse(trials, op)

    op = set_default(op, 'trace_width', 1.5);
    op = set_default(op, 'cmapname', 'jet');
    op = set_default(op, 'smooth_timecourses', true);
    op = set_default(op, 'smooth_method', 'gaussian');
    op = set_default(op, 'smooth_windowsize', 30);
    op = set_default(op, 'condval_inds_to_plot', []);
    op = set_default(op, 'y_ax_hardlims', []);
    op = set_default(op, 'epochs', []);
    op = set_default(op, 'epochs_to_label', []);
    op = set_default(op, 'epoch_colors', []);
    op = set_default(op, 'epoch_alpha', 0.12);
    op = set_default(op, 'epoch_label_fontsize', 8);
    op = set_default(op, 'show_legend', true);
    op = set_default(op, 'legend_location', 'northeast');
    op = set_default(op, 'legend_fontsize', 6);
    op = set_default(op, 'legend_box', 'on');
    op = set_default(op, 'legend_auto_avoid_overlap', true);
    op = set_default(op, 'sem_fill_color', [0.35 0.35 0.35]);
    op = set_default(op, 'sem_fill_alpha', 0.34);

    hax = gca;

    if isfield(op, 'align_stats') && isfield(op, 'resp_grpd')
        align_stats = op.align_stats;
        resp_grpd = op.resp_grpd;
    else
        [trials, align_stats, resp_grpd, op] = ...
            sort_responses_by_condition_prealigned_single_electrode(trials, op);
    end

    global_data_max = -Inf;
    global_data_min = Inf;

    collision_time_values = [];
    collision_y_lower_values = [];
    collision_y_upper_values = [];

    hold(hax, 'on')

    %% Plot SEM bands first.
    if ~isempty(resp_grpd) && height(resp_grpd) > 0
        for icond = 1:height(resp_grpd)

            if isempty(resp_grpd.resp_mean{icond}) || isempty(resp_grpd.sem{icond})
                continue
            end

            sem_lower = resp_grpd.resp_mean{icond} - resp_grpd.sem{icond};
            sem_upper = resp_grpd.resp_mean{icond} + resp_grpd.sem{icond};

            plot_indices = resp_grpd.n_good_trials{icond} > 0;

            if nnz(plot_indices) == 0
                continue
            end

            lower_to_plot = sem_lower(plot_indices);
            upper_to_plot = sem_upper(plot_indices);

            if op.smooth_timecourses
                lower_to_plot = smoothdata(lower_to_plot, 2, op.smooth_method, op.smooth_windowsize);
                upper_to_plot = smoothdata(upper_to_plot, 2, op.smooth_method, op.smooth_windowsize);
            end

            finite_sem_values = [lower_to_plot(:); upper_to_plot(:)];
            finite_sem_values = finite_sem_values(isfinite(finite_sem_values));

            if ~isempty(finite_sem_values)
                global_data_max = max(global_data_max, max(finite_sem_values));
                global_data_min = min(global_data_min, min(finite_sem_values));
            end

            times_to_plot = align_stats.times_aligned(plot_indices);

            collision_time_values = [collision_time_values; times_to_plot(:)]; %#ok<AGROW>
            collision_y_lower_values = [collision_y_lower_values; lower_to_plot(:)]; %#ok<AGROW>
            collision_y_upper_values = [collision_y_upper_values; upper_to_plot(:)]; %#ok<AGROW>

            hfill = fill(hax, ...
                [times_to_plot, fliplr(times_to_plot)], ...
                [lower_to_plot, fliplr(upper_to_plot)], ...
                op.sem_fill_color, ...
                'HandleVisibility', 'off');

            hfill.LineStyle = 'none';
            hfill.EdgeColor = op.sem_fill_color;
            hfill.FaceAlpha = op.sem_fill_alpha;
        end
    end

    %% Determine condition indices to plot.
    if isempty(op.condval_inds_to_plot)
        requested_condition_indices = 1:height(resp_grpd);
    else
        requested_condition_indices = op.condval_inds_to_plot;
    end

    valid_condition_indices = [];

    for iidx = 1:numel(requested_condition_indices)
        this_condition_index = requested_condition_indices(iidx);

        if this_condition_index <= height(resp_grpd) && ...
                ~isempty(resp_grpd.resp_mean{this_condition_index})

            valid_condition_indices(end+1) = this_condition_index; %#ok<AGROW>
        end
    end

    if isempty(valid_condition_indices)
        text(hax, 0.5, 0.5, 'No valid trials', ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center');

        trials_out = trials;
        op_out = op;
        return
    end

    %% Plot condition mean timecourses.
    timecourses_condition_by_time = cell2mat(resp_grpd.resp_mean(valid_condition_indices));
    timecourses_time_by_condition = timecourses_condition_by_time';

    if op.smooth_timecourses
        timecourses_time_by_condition = smoothdata( ...
            timecourses_time_by_condition, ...
            1, ...
            op.smooth_method, ...
            op.smooth_windowsize);
    end

    finite_timecourse_values = timecourses_time_by_condition(isfinite(timecourses_time_by_condition));

    if ~isempty(finite_timecourse_values)
        global_data_max = max(global_data_max, max(finite_timecourse_values));
        global_data_min = min(global_data_min, min(finite_timecourse_values));
    end

    h_timecourse = plot(hax, align_stats.times_aligned, timecourses_time_by_condition);

    n_conditions_to_plot = numel(valid_condition_indices);

    cmap = colormap(hax, op.cmapname);
    color_indices = round(linspace(1, size(cmap, 1), n_conditions_to_plot));

    set(h_timecourse, ...
        {'LineWidth'}, num2cell(repmat(op.trace_width, n_conditions_to_plot, 1)));

    set(h_timecourse, ...
        {'Color'}, num2cell(cmap(color_indices, :), 2));

    for icond_plot = 1:n_conditions_to_plot
        this_line_y = timecourses_time_by_condition(:, icond_plot);

        collision_time_values = [collision_time_values; align_stats.times_aligned(:)]; %#ok<AGROW>
        collision_y_lower_values = [collision_y_lower_values; this_line_y(:)]; %#ok<AGROW>
        collision_y_upper_values = [collision_y_upper_values; this_line_y(:)]; %#ok<AGROW>
    end

    %% Zero line.
    yline(hax, 0, ...
        'LineWidth', op.yline_zero_width, ...
        'Color', op.yline_zero_color, ...
        'LineStyle', op.yline_zero_style, ...
        'HandleVisibility', 'off');

    %% Initial Axis limits.
    xlim(hax, [min(align_stats.times_aligned), max(align_stats.times_aligned)])

    if ~isempty(op.y_ax_hardlims)
        ylim(hax, op.y_ax_hardlims);
    elseif isfinite(global_data_max) && isfinite(global_data_min)
        y_range = global_data_max - global_data_min;

        if y_range == 0
            y_range = 1;
        end

        y_bottom = min(global_data_min - 0.05 * y_range, -0.02 * y_range);

        if ~isempty(op.epochs)
            y_top = global_data_max + 0.25 * y_range;
        else
            y_top = global_data_max + 0.08 * y_range;
        end

        ylim(hax, [y_bottom, y_top]);
    end

    %% Axis labels.
    xlabel(hax, 'Time (sec)')

    if isfield(op, 'resp_signal')
        ylabel(hax, ['normed ', op.resp_signal, ' power'])
    else
        ylabel(hax, 'normed power')
    end

    %% Purely Analytical Legend Placement (BEFORE drawing epoch shading).
    if op.show_legend && ~is_empty_sort_condition(op.sort_cond) && height(resp_grpd) > 0

        legend_strings = cellstr(safe_string_column(resp_grpd.condval(valid_condition_indices)));
        legend_title_string = option_value_to_title_string(op.sort_cond);

        hleg = legend(hax, h_timecourse, legend_strings, ...
            'Interpreter', 'none', ...
            'Location', op.legend_location, ...
            'FontSize', op.legend_fontsize, ...
            'Box', op.legend_box);

        title(hleg, legend_title_string, 'Interpreter', 'none');

        try
            hleg.AutoUpdate = 'off';
        catch
        end

        if op.legend_auto_avoid_overlap
            place_legend_analytically( ...
                hax, ...
                hleg, ...
                op, ...
                collision_time_values, ...
                collision_y_lower_values, ...
                collision_y_upper_values);
        end
    end

    %% Epoch backgrounds and labels (PLOTTED STRICTLY ONCE AFTER YLIM IS FINALIZED).
    if ~isempty(op.epochs)
        plot_epochs_on_timecourse(hax, op);
    end

    box(hax, 'off')
    hold(hax, 'off')

    trials_out = trials;
    op_out = op;
end

%% ========================================================================
% Helper: analytically place legend in matrix space (instantaneous)
% ========================================================================

function place_legend_analytically(hax, hleg, op, data_x, data_y_low, data_y_high)

    if isempty(data_x) || isempty(data_y_low) || isempty(data_y_high)
        return
    end

    % Get approximate legend width and height in normalized coordinates
    leg_units = hleg.Units;
    hleg.Units = 'normalized';
    leg_pos = hleg.Position;
    hleg.Units = leg_units;

    w_leg = max(leg_pos(3), 0.22);
    h_leg = max(leg_pos(4), 0.16);
    margin = 0.02;

    % Define standard inside corners: [x_min, y_min, width, height]
    candidate_locations = {'northeast', 'northwest', 'southeast', 'southwest'};
    candidate_boxes = [ ...
        1 - margin - w_leg,  1 - margin - h_leg,  w_leg, h_leg;  % northeast
        margin,              1 - margin - h_leg,  w_leg, h_leg;  % northwest
        1 - margin - w_leg,  margin,              w_leg, h_leg;  % southeast
        margin,              margin,              w_leg, h_leg;  % southwest
    ];

    x_limits = xlim(hax);
    y_limits = ylim(hax);

    finite_mask = isfinite(data_x) & isfinite(data_y_low) & isfinite(data_y_high);
    dx = data_x(finite_mask);
    dy_low = data_y_low(finite_mask);
    dy_high = data_y_high(finite_mask);

    if isempty(dx)
        return
    end

    % Convert data values into normalized [0, 1] axes coordinates
    x_norm = (dx - x_limits(1)) / diff(x_limits);
    y_low_norm = (dy_low - y_limits(1)) / diff(y_limits);
    y_high_norm = (dy_high - y_limits(1)) / diff(y_limits);

    % Precompute normalized centers of epoch text labels
    epoch_label_x_norm = [];

    if ~isempty(op.epochs)
        epoch_durations = op.epochs.dur_fix(:);
        epoch_starts = [0; cumsum(epoch_durations(1:end-1))];
        epoch_mids = epoch_starts + epoch_durations / 2;
        epoch_label_x_norm = (epoch_mids - x_limits(1)) / diff(x_limits);
    end

    % Check candidates purely in RAM
    for ic = 1:size(candidate_boxes, 1)

        box = candidate_boxes(ic, :);
        box_x1 = box(1); box_x2 = box(1) + box(3);
        box_y1 = box(2); box_y2 = box(2) + box(4);

        % 1. Data collision check
        in_x = (x_norm >= box_x1 & x_norm <= box_x2);
        data_collide = any(in_x & (y_high_norm >= box_y1 & y_low_norm <= box_y2));

        if data_collide
            continue
        end

        % 2. Epoch label collision check (labels reside in top region y_norm >= 0.82)
        label_collide = false;

        if ~isempty(epoch_label_x_norm) && box_y2 >= 0.82
            for el_x = epoch_label_x_norm'
                if (el_x + 0.05 >= box_x1) && (el_x - 0.05 <= box_x2)
                    label_collide = true;
                    break
                end
            end
        end

        if ~label_collide
            hleg.Location = candidate_locations{ic};
            return
        end
    end

    % If all 4 corners collide with data/labels, compute required upper y-limit expansion
    max_data_y = max(dy_high);
    available_top_frac = 1.0 - h_leg - 0.06;

    if available_top_frac > 0.1
        new_y_max = y_limits(1) + (max_data_y - y_limits(1)) / available_top_frac;
        ylim(hax, [y_limits(1), new_y_max]);
    end

    hleg.Location = 'northeast';
end

function candidate_locations = normalize_legend_candidate_locations(candidate_locations)

    if isstring(candidate_locations)
        candidate_locations = cellstr(candidate_locations);
    elseif ischar(candidate_locations)
        candidate_locations = {candidate_locations};
    elseif iscell(candidate_locations)
        candidate_locations = cellstr(safe_string_column(candidate_locations));
    else
        error('op.legend_candidate_locations must be a char, string array, or cell array.')
    end
end

function output_cell = unique_stable_cellstr(input_cell)

    input_strings = string(input_cell);
    [~, unique_indices] = unique(input_strings, 'stable');
    output_cell = cellstr(input_strings(sort(unique_indices)));
end

%% ========================================================================
% Helper: plot epoch backgrounds (Executes strictly ONCE per subplot)
% ========================================================================

function plot_epochs_on_timecourse(hax, op)

    % Clear any pre-existing epoch elements to prevent alpha stack-up
    delete(findobj(hax, 'Tag', 'epoch_graphic'));

    epochs = op.epochs;

    if isempty(op.epoch_colors)
        op.epoch_colors = lines(height(epochs));
    end

    if size(op.epoch_colors, 1) < height(epochs)
        error('op.epoch_colors must have at least one RGB row per epoch.')
    end

    if ~isempty(op.epochs_to_label)
        if ~iscell(op.epochs_to_label) && ~isstring(op.epochs_to_label)
            error('op.epochs_to_label must be empty, a cell array, or a string array.')
        end

        epochs_to_show = cellstr(safe_string_column(op.epochs_to_label));
    else
        epochs_to_show = cellstr(safe_string_column(epochs.epoch));
    end

    epoch_durations = epochs.dur_fix(:);

    epoch_starts = [0; cumsum(epoch_durations(1:end-1))];
    epoch_ends = epoch_starts + epoch_durations;
    epoch_mids = epoch_starts + epoch_durations / 2;

    ylims = ylim(hax);
    y_range = ylims(2) - ylims(1);

    hold(hax, 'on')

    for iepoch = 1:height(epochs)
        x_patch = [epoch_starts(iepoch), epoch_ends(iepoch), ...
                   epoch_ends(iepoch), epoch_starts(iepoch)];

        y_patch = [ylims(1), ylims(1), ylims(2), ylims(2)];

        h_patch = patch(hax, ...
            x_patch, ...
            y_patch, ...
            op.epoch_colors(iepoch, :), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', op.epoch_alpha, ...
            'HandleVisibility', 'off', ...
            'Tag', 'epoch_graphic');

        uistack(h_patch, 'bottom');
    end

    for iepoch = 1:height(epochs)
        xline(hax, epoch_starts(iepoch), ...
            'LineStyle', op.epoch_divider_linestyle, ...
            'LineWidth', op.epoch_divider_linewidth, ...
            'Color', [0.5 0.5 0.5], ...
            'HandleVisibility', 'off', ...
            'Tag', 'epoch_graphic');
    end

    xline(hax, epoch_ends(end), ...
        'LineStyle', op.epoch_divider_linestyle, ...
        'LineWidth', op.epoch_divider_linewidth, ...
        'Color', [0.5 0.5 0.5], ...
        'HandleVisibility', 'off', ...
        'Tag', 'epoch_graphic');

    epoch_names = safe_string_column(epochs.epoch);

    for ilabel = 1:numel(epochs_to_show)

        epoch_name = string(epochs_to_show{ilabel});
        epoch_index = find(epoch_names == epoch_name, 1);

        if isempty(epoch_index)
            continue
        end

        x_text = epoch_mids(epoch_index);
        y_text = ylims(2) - 0.03 * y_range;

        text(hax, x_text, y_text, epoch_name, ...
            'FontSize', op.epoch_label_fontsize, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'none', ...
            'Margin', 1, ...
            'Tag', 'epoch_graphic');
    end
end

%% ========================================================================
% General utilities
% ========================================================================

function op = set_default(op, field_name, default_value)

    if ~isfield(op, field_name) || isempty(op.(field_name))
        op.(field_name) = default_value;
    end
end

function tf = is_empty_sort_condition(sort_condition)

    if isempty(sort_condition)
        tf = true;
        return
    end

    if iscell(sort_condition)
        tf = isempty(sort_condition{1}) || strcmp(char(string(sort_condition{1})), '');
        return
    end

    tf = strcmp(char(string(sort_condition)), '');
end

function response_matrix = normalize_electrode_timecourse(raw_timecourse, expected_number_of_trials)

    if iscell(raw_timecourse)
        if isrow(raw_timecourse)
            raw_timecourse = raw_timecourse';
        end

        response_matrix = cell2mat(raw_timecourse);
    else
        response_matrix = raw_timecourse;
    end

    if isvector(response_matrix)
        response_matrix = response_matrix(:)';
    end

    if size(response_matrix, 1) ~= expected_number_of_trials && ...
            size(response_matrix, 2) == expected_number_of_trials

        response_matrix = response_matrix';
    end
end

function time_vector = extract_resp_time_vector(resp, row_index)

    time_vector = [];

    if ~ismember('times', resp.Properties.VariableNames)
        return
    end

    raw_times = resp.times;

    if iscell(raw_times)
        candidate = raw_times{row_index};
    else
        candidate = raw_times(row_index, :);
    end

    if isempty(candidate)
        return
    end

    if iscell(candidate)
        if isempty(candidate)
            return
        end
        candidate = candidate{1};
    end

    if isnumeric(candidate)
        time_vector = candidate(:)';
    end
end

function time_vector = extract_trial_time_vector(trials)

    time_vector = [];

    if ~ismember('times', trials.Properties.VariableNames)
        return
    end

    if iscell(trials.times)
        if ~isempty(trials.times{1})
            time_vector = trials.times{1};
        end
    elseif isnumeric(trials.times)
        time_vector = trials.times(1, :);
    end

    if ~isempty(time_vector)
        time_vector = time_vector(:)';
    end
end

function resp_grpd = make_single_condition_response_group(condition_name, resp_aligned, align_stats)

    empty_cell_column = cell(1, 1);

    resp_grpd = table({condition_name}, ...
                      empty_cell_column, ...
                      empty_cell_column, ...
                      'VariableNames', {'condval', 'resp', 'resp_mean'});

    resp_grpd.resp{1} = resp_aligned;
    resp_grpd.resp_mean{1} = align_stats.mean;
    resp_grpd.std{1} = align_stats.std;
    resp_grpd.n_good_trials{1} = align_stats.n_good_trials;
    resp_grpd.sem{1} = align_stats.sem;
end

function region_definition_table = make_region_definition_table(resp, requested_region_order)

    if ~ismember('region', resp.Properties.VariableNames)
        region_definition_table = table();
        return
    end

    response_region_names = safe_string_column(resp.region);

    if ~isempty(requested_region_order)
        candidate_regions = safe_string_column(requested_region_order);
    else
        candidate_regions = unique(response_region_names, 'stable');
    end

    candidate_regions = candidate_regions(candidate_regions ~= "");

    region_names = {};
    region_counts = [];

    for iregion = 1:numel(candidate_regions)
        this_region = candidate_regions(iregion);
        this_count = nnz(response_region_names == this_region);

        if this_count > 0
            region_names{end+1, 1} = char(this_region); %#ok<AGROW>
            region_counts(end+1, 1) = this_count; %#ok<AGROW>
        end
    end

    region_definition_table = table(region_names, region_counts, ...
        'VariableNames', {'region', 'n_elcs'});
end

function hfig = create_positioned_figure(figure_index, op)

    monitor_positions = get(groot, 'MonitorPositions');
    n_monitors = size(monitor_positions, 1);

    if n_monitors == 0
        hfig = figure('Color', 'w');
        return
    end

    if isfield(op, 'screen_order') && ~isempty(op.screen_order)
        screen_order = op.screen_order(:)';

        if ~isnumeric(screen_order) || any(screen_order ~= round(screen_order))
            error('op.screen_order must be a numeric vector of integer screen IDs.')
        end

        if any(screen_order < 1) || any(screen_order > n_monitors)
            error('op.screen_order contains a screen ID outside the valid range 1:%d.', n_monitors)
        end
    else
        if n_monitors == 1
            screen_order = 1;
        else
            screen_order = [n_monitors, 1:(n_monitors-1)];
        end
    end

    screen_order_index = mod(figure_index - 1, numel(screen_order)) + 1;
    monitor_id = screen_order(screen_order_index);

    monitor_position = monitor_positions(monitor_id, :);

    hfig = figure( ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'OuterPosition', monitor_position);

    drawnow

    try
        hfig.WindowState = 'maximized';
    catch
    end
end

function title_string = make_electrode_title(resp, row_index)

    if ismember('sub', resp.Properties.VariableNames)
        all_subject_strings = safe_string_column(resp.sub);
        subject_string = char(all_subject_strings(row_index));

        if isempty(subject_string)
            subject_string = 'unknown-sub';
        end
    else
        subject_string = 'unknown-sub';
    end

    if ismember('chan', resp.Properties.VariableNames)
        channel_value = resp.chan(row_index);

        if iscell(channel_value)
            channel_string = char(safe_single_value_to_string(channel_value{1}));
        else
            channel_string = char(safe_single_value_to_string(channel_value));
        end

        if isempty(channel_string)
            channel_string = 'unknown-chan';
        end
    else
        channel_string = 'unknown-chan';
    end

    if ismember('region', resp.Properties.VariableNames)
        all_region_strings = safe_string_column(resp.region);
        region_string = char(all_region_strings(row_index));

        if isempty(region_string)
            region_string = 'unknown-region';
        end
    else
        region_string = 'unknown-region';
    end

    title_string = sprintf('%s | chan %s | %s', subject_string, channel_string, region_string);
end

function title_string = option_value_to_title_string(option_value)

    if iscell(option_value)
        if numel(option_value) >= 2 && isnumeric(option_value{2})
            title_string = sprintf('%s{%d}', char(string(option_value{1})), option_value{2});
        else
            title_string = strjoin(cellstr(safe_string_column(option_value)), ', ');
        end
    else
        title_string = char(safe_single_value_to_string(option_value));
    end
end

function string_column = safe_string_column(input_values)

    if isstring(input_values)
        string_column = input_values(:);
        string_column(ismissing(string_column)) = "";
        return
    end

    if ischar(input_values)
        string_column = string(cellstr(input_values));
        string_column = string_column(:);
        string_column(ismissing(string_column)) = "";
        return
    end

    if iscategorical(input_values)
        string_column = string(input_values(:));
        string_column(ismissing(string_column)) = "";
        return
    end

    if isnumeric(input_values) || islogical(input_values)
        string_column = string(input_values(:));
        string_column(ismissing(string_column)) = "";
        return
    end

    if iscell(input_values)
        string_column = strings(numel(input_values), 1);

        for ivalue = 1:numel(input_values)
            string_column(ivalue) = safe_single_value_to_string(input_values{ivalue});
        end

        string_column(ismissing(string_column)) = "";
        return
    end

    try
        string_column = string(input_values(:));
        string_column(ismissing(string_column)) = "";
    catch
        error('Could not safely convert variable of class "%s" to string.', class(input_values))
    end
end

function output_string = safe_single_value_to_string(input_value)

    if isempty(input_value)
        output_string = "";
        return
    end

    if isstring(input_value)
        input_value = input_value(:);
        input_value(ismissing(input_value)) = [];

        if isempty(input_value)
            output_string = "";
        elseif numel(input_value) == 1
            output_string = input_value;
        else
            output_string = strjoin(input_value, "|");
        end

        return
    end

    if ischar(input_value)
        output_string = string(input_value);
        return
    end

    if iscategorical(input_value)
        output_string = safe_single_value_to_string(string(input_value));
        return
    end

    if isnumeric(input_value) || islogical(input_value)
        if isempty(input_value)
            output_string = "";
        elseif isscalar(input_value)
            output_string = string(input_value);
        else
            output_string = strjoin(string(input_value(:)), "|");
        end

        return
    end

    if iscell(input_value)
        if isempty(input_value)
            output_string = "";
            return
        end

        nested_strings = strings(numel(input_value), 1);

        for jvalue = 1:numel(input_value)
            nested_strings(jvalue) = safe_single_value_to_string(input_value{jvalue});
        end

        nested_strings = nested_strings(nested_strings ~= "");

        if isempty(nested_strings)
            output_string = "";
        elseif numel(nested_strings) == 1
            output_string = nested_strings;
        else
            output_string = strjoin(nested_strings, "|");
        end

        return
    end

    try
        converted_value = string(input_value);
        converted_value = converted_value(:);
        converted_value(ismissing(converted_value)) = [];

        if isempty(converted_value)
            output_string = "";
        elseif numel(converted_value) == 1
            output_string = converted_value;
        else
            output_string = strjoin(converted_value, "|");
        end
    catch
        output_string = "";
    end
end