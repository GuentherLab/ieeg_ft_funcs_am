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
    field_default('op', 'leg_pos_adjust', 0.22)
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
    elseif strcmp(op.sort_cond, '')
        % No sorting - create minimal align_stats
        [trials, align_stats, op] = minimal_align_prealigned(trials, op);
        resp_grpd = table();
    end
    
    %% Plot error bands (SEM)
    if ~strcmp(op.sort_cond, '') && height(resp_grpd) > 0
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
        if ~strcmp(op.sort_cond, '') && height(resp_grpd) > 0
            op.condval_inds_to_plot = 1:height(resp_grpd);
        else
            op.condval_inds_to_plot = 1;
        end
    end
    
    nvals_to_plot = length(op.condval_inds_to_plot);
    
    if ~strcmp(op.sort_cond, '') && height(resp_grpd) > 0
        timecourses_to_plot = cell2mat(resp_grpd.resp_mean(op.condval_inds_to_plot, :))';
        legend_strs = resp_grpd.condval(op.condval_inds_to_plot);
    else
        timecourses_to_plot = nanmean(trials.resp_aligned, 1);
        legend_strs = {'All Trials'};
    end
    
    if op.smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 1, op.smooth_method, op.smooth_windowsize);
    end
    
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
    
    %% Plot epoch backgrounds and labels
    if ~isempty(op.epochs)
        plot_epochs_on_timecourse(hax, op);
    end
    
    %% Formatting
    xlim(hax, [min(align_stats.times_aligned), max(align_stats.times_aligned)])
    if ~isempty(op.y_ax_hardlims)
        ylim(hax, op.y_ax_hardlims)
    end
    
    xlabel(hax, 'Time (sec)')
    if isfield(op, 'resp_signal')
        ylabel(hax, ['normed ', op.resp_signal, ' power'])
    else
        ylabel(hax, 'normed power')
    end
    
    set(gcf, 'Color', [1 1 1])
    
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
function [trials_out, align_stats, op_out] = minimal_align_prealigned(trials, op)
    % For pre-warped data with no condition sorting
    
    field_default('op', 'samp_period', nan)
    
    if ismember('resp_aligned', trials.Properties.VariableNames)
        resp_aligned = trials.resp_aligned;
    else
        resp_aligned = cell2mat(trials.resp_unaligned);
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
    
    % Convert to trial-relative (start at 0)
    times_aligned = times_aligned - times_aligned(1);
    
    align_stats.times_aligned = times_aligned;
    align_stats.n_tpoints_pre_fixed = length(times_aligned);
    align_stats.n_tpoints_post_fixed = 0;
    align_stats.samp_period = op.samp_period;
    align_stats.mean = nanmean(resp_aligned, 1);
    align_stats.std = nanstd(resp_aligned);
    align_stats.sem = align_stats.std ./ sqrt(sum(~isnan(resp_aligned), 1));
    
    trials_out = trials;
    op_out = op;
end

%% ========== HELPER: Plot epochs as backgrounds with labels and dividers ==========
function plot_epochs_on_timecourse(hax, op)
    % Plot epoch backgrounds, dividers, and text labels using axis coordinates only
    
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
    
    %% Get current axis limits
    hold(hax, 'on')
    ylims = ylim(hax);
    xlims = xlim(hax);
    y_range = ylims(2) - ylims(1);
    
    %% Draw epoch backgrounds using patch
    for iepoch = 1:height(epochs)
        % Create patch for epoch background
        x_patch = [epoch_starts(iepoch), epoch_ends(iepoch), ...
                   epoch_ends(iepoch), epoch_starts(iepoch)];
        y_patch = [ylims(1), ylims(1), ylims(2), ylims(2)];
        
        h_patch = patch(hax, x_patch, y_patch, op.epoch_colors(iepoch, :), ...
            'EdgeColor', 'none', 'FaceAlpha', op.epoch_alpha);
        set(h_patch, 'HandleVisibility', 'off');
        uistack(h_patch, 'bottom');  % Send to back
    end
    
    %% Draw epoch dividers (vertical lines at epoch boundaries)
    for iepoch = 1:height(epochs)
        xline(hax, epoch_starts(iepoch), 'LineStyle', op.epoch_divider_linestyle, ...
            'LineWidth', op.epoch_divider_linewidth, 'Color', [0.5 0.5 0.5], ...
            'HandleVisibility', 'off');
    end
    % Draw final epoch end boundary
    xline(hax, epoch_ends(end), 'LineStyle', op.epoch_divider_linestyle, ...
        'LineWidth', op.epoch_divider_linewidth, 'Color', [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    
    %% Add text labels for selected epochs - INSIDE PLOT AREA
    % Position text labels at the top of each epoch, inside the plot
    for idx_label = 1:length(epochs_to_show)
        epoch_name = epochs_to_show{idx_label};
        
        % Find index in epochs table
        epoch_idx = find(strcmp(epochs.epoch, epoch_name));
        if isempty(epoch_idx)
            continue
        end
        
        % Position text at epoch midpoint, near top of plot
        x_text = epoch_mids(epoch_idx);
        y_text = ylims(2) - 0.08 * y_range;  % 92% up from bottom
        
        % Create text object (no annotation box needed)
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
    
    % Convert to trial-relative (start at 0)
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
        trials.sort_cond = trials{:, sort_col};
    end
    
    %% Handle sort_cond_vals: filter and order
    if ~isempty(op.sort_cond_vals)
        % Convert sort_cond_vals and table conditions to string for robust comparison
        cond_vals_str = cellstr(string(op.sort_cond_vals));
        
        % Filter trials and obtain exact index mapping matching op.sort_cond_vals order
        [keep_trials, trial_cond_ind] = ismember(string(trials.sort_cond), string(cond_vals_str));
        
        % Apply filter to trials and response matrix
        trials = trials(keep_trials, :);
        resp_aligned = resp_aligned(keep_trials, :);
        trial_cond_ind = trial_cond_ind(keep_trials);
        
        op.sort_cond_vals = cond_vals_str;
    else
        % If not specified, get all unique values
        unique_vals = unique(trials.sort_cond);
        if isnumeric(unique_vals)
            unique_vals = unique_vals(~isnan(unique_vals));
        end
        op.sort_cond_vals = cellstr(string(unique_vals));
        
        [~, trial_cond_ind] = ismember(string(trials.sort_cond), string(op.sort_cond_vals));
    end
    
    nconds = length(op.sort_cond_vals);
    
    %% Group responses by condition
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
            % Empty condition
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