function [trials_out, resp_grpd, align_stats, op_out] = plot_resp_timecourse(trials, op)
    %% PLOT_RESP_TIMECOURSE - Plot electrode timecourses with optional epoch backgrounds
    % 
    % For use with PRE-WARPED data where all trials have identical timepoints
    %
    % Inputs:
    %   trials: table with one row per trial, including sort_cond variable
    %   op: options struct
    %
    % Outputs:
    %   trials_out, resp_grpd, align_stats, op_out
    
    %% Set defaults
    field_default('op', 'sort_cond', '')
    field_default('op', 'do_condition_sorting', 1)
    field_default('op', 'trace_width', 1.5)
    field_default('op', 'cmapname', 'jet')
    field_default('op', 'smooth_timecourses', 1)
    field_default('op', 'smooth_method', 'gaussian')
    field_default('op', 'smooth_windowsize', 30)
    field_default('op', 'condval_inds_to_plot', [])
    field_default('op', 'newfig', 1)
    field_default('op', 'yline_zero_width', 0.25)
    field_default('op', 'yline_zero_color', [0.8 0.8 0.8])
    field_default('op', 'yline_zero_style', '-')
    field_default('op', 'y_ax_hardlims', [])
    field_default('op', 'epochs', [])
    field_default('op', 'epochs_to_label', [])
    field_default('op', 'epoch_colors', [])
    field_default('op', 'epoch_alpha', 0.12)
    field_default('op', 'epoch_label_height', 0.85)
    field_default('op', 'epoch_label_fontsize', 7)
    field_default('op', 'epoch_divider_linestyle', ':')
    field_default('op', 'epoch_divider_linewidth', 1.5)
    field_default('op', 'samp_period', nan)
    field_default('op', 'sort_cond_vals', {})
    
    %% Create figure if needed
    if op.newfig
        hfig = figure('Color', [1 1 1]);
        box off
    end
    hax = gca;
    
    %% Handle sorting
    if op.do_condition_sorting && ~strcmp(op.sort_cond, '')
        if iscell(op.sort_cond)
            trials.sort_cond = trials{:, op.sort_cond{1}}(:, op.sort_cond{2});
        else
            trials.sort_cond = trials{:, op.sort_cond};
        end
        [trials, align_stats, resp_grpd, op] = sort_responses_by_condition_prealigned(trials, op);
    elseif ~op.do_condition_sorting
        % Pre-computed grouping - use provided align_stats
        align_stats = op.align_stats;
        resp_grpd = op.resp_grpd;
    else
        % No sorting or empty sort_cond - generate overall mean and SEM across all trials
        [trials, align_stats, resp_grpd, op] = minimal_align_prealigned(trials, op);
    end
    
    %% Track peak and trough across all signals + SEM bounds
    global_data_max = -Inf;
    global_data_min = Inf;
    
    %% Plot error bands (SEM)
    if ~isempty(resp_grpd) && height(resp_grpd) > 0
        nconds = height(resp_grpd);
        for icond = 1:nconds
            this_cond_sem_lims = [resp_grpd.resp_mean{icond} - resp_grpd.sem{icond};
                                  resp_grpd.resp_mean{icond} + resp_grpd.sem{icond}];
            plotinds = resp_grpd.n_good_trials{icond} > 0;
            
            if nnz(plotinds) > 0
                lowlims = this_cond_sem_lims(1, plotinds);
                uplims = fliplr(this_cond_sem_lims(2, plotinds));
                
                if op.smooth_timecourses
                    lowlims = smoothdata(lowlims, 2, op.smooth_method, op.smooth_windowsize);
                    uplims = smoothdata(uplims, 2, op.smooth_method, op.smooth_windowsize);
                end
                
                global_data_max = max([global_data_max, max(uplims), max(lowlims)]);
                global_data_min = min([global_data_min, min(uplims), min(lowlims)]);
                
                times_plot = align_stats.times_aligned(plotinds);
                hfill = fill(hax, [times_plot, fliplr(times_plot)], [lowlims, uplims], ...
                    [0.8 0.8 0.8], 'HandleVisibility', 'off');
                hfill.LineStyle = 'none';
                hfill.EdgeColor = [0.8 0.8 0.8];
                hold(hax, 'on')
            end
        end
    end
    
    %% Plot timecourses
    if isempty(op.condval_inds_to_plot)
        if ~isempty(resp_grpd) && height(resp_grpd) > 0
            op.condval_inds_to_plot = 1:height(resp_grpd);
        else
            op.condval_inds_to_plot = 1;
        end
    end
    
    nvals_to_plot = length(op.condval_inds_to_plot);
    
    if ~isempty(resp_grpd) && height(resp_grpd) > 0
        timecourses_to_plot = cell2mat(resp_grpd.resp_mean(op.condval_inds_to_plot, :))';
        legend_strs = resp_grpd.condval(op.condval_inds_to_plot);
    else
        timecourses_to_plot = nanmean(trials.resp_aligned, 1);
        legend_strs = {'All Trials'};
    end
    
    if op.smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 1, op.smooth_method, op.smooth_windowsize);
    end
    
    global_data_max = max([global_data_max, max(timecourses_to_plot(:))]);
    global_data_min = min([global_data_min, min(timecourses_to_plot(:))]);
    
    hold(hax, 'on')
    h_timecourse = plot(hax, align_stats.times_aligned, timecourses_to_plot);
    
    % Apply colormap
    cmap = colormap(hax, op.cmapname);
    color_indices = round(linspace(1, size(cmap, 1), nvals_to_plot));
    set(h_timecourse, {'LineWidth'}, num2cell(repmat(op.trace_width, nvals_to_plot, 1)));
    set(h_timecourse, {'Color'}, num2cell(cmap(color_indices, :), 2));
    
    %% Add zero line
    yline(hax, 0, 'LineWidth', op.yline_zero_width, 'Color', op.yline_zero_color, ...
        'LineStyle', op.yline_zero_style);
    
    %% Dynamic Y-Axis Limits Adjustment (Headroom for Epoch Labels)
    xlim(hax, [min(align_stats.times_aligned), max(align_stats.times_aligned)])
    
    if ~isempty(op.y_ax_hardlims)
        ylim(hax, op.y_ax_hardlims);
    elseif isfinite(global_data_max) && isfinite(global_data_min)
        y_range = global_data_max - global_data_min;
        if y_range == 0, y_range = 1; end
        
        % Set bottom limit (preserve baseline comfortably)
        y_bottom = min(global_data_min - 0.05 * y_range, -0.02 * y_range);
        
        % If epochs are labeled, add headroom above highest peak/SEM for text labels
        if ~isempty(op.epochs)
            headroom_factor = 0.28; 
            y_top = global_data_max + y_range * headroom_factor;
        else
            y_top = global_data_max + 0.08 * y_range;
        end
        
        ylim(hax, [y_bottom, y_top]);
    end
    
    %% Plot epoch backgrounds and labels AFTER y-limits are finalized
    if ~isempty(op.epochs)
        plot_epochs_on_timecourse(hax, op);
    end
    
    %% Formatting
    xlabel(hax, 'Time (sec)')
    if isfield(op, 'resp_signal')
        ylabel(hax, ['normed ', op.resp_signal, ' power'])
    else
        ylabel(hax, 'normed power')
    end
    
    set(gcf, 'Color', [1 1 1])
    
    % Only show legend when actively sorting by condition
    if ~strcmp(op.sort_cond, '') && height(resp_grpd) > 0
        hleg = legend(hax, legend_strs{:}, 'Interpreter', 'none', 'Location', 'bestoutside');
        title(hleg, op.sort_cond)
    end
    
    hold(hax, 'off')
    box(hax, 'off')
    
    %% Output
    trials_out = trials;
    op_out = op;
end

%% ========== HELPER: Minimal alignment for pre-warped data ==========
function [trials_out, align_stats, resp_grpd, op_out] = minimal_align_prealigned(trials, op)
    % For pre-warped data with no condition sorting (or empty sort_cond)
    
    field_default('op', 'samp_period', nan)
    
    if ismember('resp_aligned', trials.Properties.VariableNames)
        resp_aligned = trials.resp_aligned;
    else
        resp_aligned = cell2mat(trials.resp_unaligned);
    end
    
    if iscell(resp_aligned)
        resp_aligned = cell2mat(resp_aligned);
    end
    
    if size(resp_aligned, 1) == 1
        resp_aligned = resp_aligned';
    end
    
    if isnan(op.samp_period)
        if iscell(trials.times)
            first_times = trials.times{1};
        else
            first_times = trials.times(1, :);
        end
        op.samp_period = mean(diff(first_times));
    end
    
    if iscell(trials.times)
        times_aligned = trials.times{1};
    else
        times_aligned = trials.times(1, :);
    end
    
    times_aligned = times_aligned - times_aligned(1);
    
    align_stats.times_aligned = times_aligned;
    align_stats.n_tpoints_pre_fixed = length(times_aligned);
    align_stats.n_tpoints_post_fixed = 0;
    align_stats.samp_period = op.samp_period;
    align_stats.mean = nanmean(resp_aligned, 1);
    align_stats.std = nanstd(resp_aligned);
    align_stats.sem = align_stats.std ./ sqrt(sum(~isnan(resp_aligned), 1));
    align_stats.sem_lims = [align_stats.mean + align_stats.sem; 
                            align_stats.mean - align_stats.sem];
    
    % Populate 1-row resp_grpd for 'All Trials' so error bars and means plot seamlessly
    celcol = cell(1, 1);
    resp_grpd = table({'All Trials'}, celcol, celcol, ...
        'VariableNames', {'condval', 'resp', 'resp_mean'});
    resp_grpd.resp{1} = resp_aligned;
    resp_grpd.resp_mean{1} = align_stats.mean;
    resp_grpd.std{1} = align_stats.std;
    resp_grpd.n_good_trials{1} = sum(~isnan(resp_aligned), 1);
    resp_grpd.sem{1} = align_stats.sem;
    
    trials_out = trials;
    op_out = op;
end

%% ========== HELPER: Plot epochs as backgrounds with labels and dividers ==========
function plot_epochs_on_timecourse(hax, op)
    epochs = op.epochs;
    
    %% Validate epochs_to_label
    if ~isempty(op.epochs_to_label)
        if ~iscell(op.epochs_to_label)
            error('op.epochs_to_label must be empty or a cell array of epoch names')
        end
        
        available_labels = epochs.epoch;
        for i = 1:length(op.epochs_to_label)
            if ~ismember(op.epochs_to_label{i}, available_labels)
                error('Epoch "%s" in op.epochs_to_label not found in epochs table', op.epochs_to_label{i})
            end
        end
        epochs_to_show = op.epochs_to_label;
    else
        epochs_to_show = epochs.epoch;
    end
    
    %% Build epoch timing (trial-relative, starting at 0)
    epoch_starts = [0; cumsum(epochs.dur_fix(1:end-1))];
    epoch_ends = epoch_starts + epochs.dur_fix;
    epoch_mids = epoch_starts + epochs.dur_fix/2;
    
    %% Get current axis limits (which now include headroom)
    hold(hax, 'on')
    ylims = ylim(hax);
    y_range = ylims(2) - ylims(1);
    
    %% Draw epoch backgrounds using patch
    for iepoch = 1:height(epochs)
        x_patch = [epoch_starts(iepoch), epoch_ends(iepoch), ...
                   epoch_ends(iepoch), epoch_starts(iepoch)];
        y_patch = [ylims(1), ylims(1), ylims(2), ylims(2)];
        
        h_patch = patch(hax, x_patch, y_patch, op.epoch_colors(iepoch, :), ...
            'EdgeColor', 'none', 'FaceAlpha', op.epoch_alpha);
        set(h_patch, 'HandleVisibility', 'off');
        uistack(h_patch, 'bottom');  % Send to back behind timecourse lines
    end
    
    %% Draw epoch dividers (vertical lines at epoch boundaries)
    for iepoch = 1:height(epochs)
        xline(hax, epoch_starts(iepoch), 'LineStyle', op.epoch_divider_linestyle, ...
            'LineWidth', op.epoch_divider_linewidth, 'Color', [0.5 0.5 0.5], ...
            'HandleVisibility', 'off');
    end
    xline(hax, epoch_ends(end), 'LineStyle', op.epoch_divider_linestyle, ...
        'LineWidth', op.epoch_divider_linewidth, 'Color', [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    
    %% Add text labels for selected epochs inside top headroom margin
    for idx_label = 1:length(epochs_to_show)
        epoch_name = epochs_to_show{idx_label};
        
        epoch_idx = find(strcmp(epochs.epoch, epoch_name));
        if isempty(epoch_idx)
            continue
        end
        
        x_text = epoch_mids(epoch_idx);
        y_text = ylims(2) - 0.03 * y_range;  
        
        text(hax, x_text, y_text, epoch_name, ...
            'FontSize', op.epoch_label_fontsize, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'none', ...
            'Margin', 1);
    end
    
    hold(hax, 'off')
end

%% ========== HELPER: Sort by condition (pre-warped) ==========
function [trials_out, align_stats, resp_grpd, op_out] = ...
    sort_responses_by_condition_prealigned(trials, op)
    
    field_default('op', 'sort_cond_vals', {})
    field_default('op', 'samp_period', nan)
    field_default('op', 'sort_cond', '')
    
    resp_aligned = trials.resp_aligned;
    
    if iscell(resp_aligned)
        resp_aligned = cell2mat(resp_aligned);
    end
    
    if size(resp_aligned, 1) == 1
        resp_aligned = resp_aligned';
    end
    
    if ismember('times', trials.Properties.VariableNames)
        if iscell(trials.times)
            times_aligned = trials.times{1};
        else
            times_aligned = trials.times(1, :);
        end
    else
        times_aligned = 1:size(resp_aligned, 2);
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
    
    if strcmp(op.sort_cond, '')
        % If sort_cond is empty, return overall statistics as a single group
        celcol = cell(1, 1);
        resp_grpd = table({'All Trials'}, celcol, celcol, ...
            'VariableNames', {'condval', 'resp', 'resp_mean'});
        resp_grpd.resp{1} = resp_aligned;
        resp_grpd.resp_mean{1} = align_stats.mean;
        resp_grpd.std{1} = align_stats.std;
        resp_grpd.n_good_trials{1} = sum(~isnan(resp_aligned), 1);
        resp_grpd.sem{1} = align_stats.sem;
        
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
        [keep_trials, trial_cond_ind] = ismember(string(trials.sort_cond), string(cond_vals_str));
        
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