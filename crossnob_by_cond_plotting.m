function op_out = crossnob_by_cond_plotting(resp, subs, op)
% CROSSNOB_BY_COND_PLOTTING Plots Crossnobis distance/SNR results by ROI and subject.
%
% Inputs:
%   resp - Response table containing 'cnob', 'sub', 'chan', 'region', etc.
%   subs - Table containing a column/variable named 'sub' specifying subjects to plot.
%   op   - Struct with optional fields:
%          .var_to_plot                  - 'mean_dist' (default) or 'snr'
%          .sort_cond_vals               - Cellstr of conditions to plot (default: all)
%          .epochs_to_analyze            - Cellstr of epochs to plot (default: all)
%          .analyze_responsive_elcs_only - (boolean) keep only responsive elcs (default: true)
%          .analyze_tuned_elcs_only      - (boolean) keep only tuned elcs (default: false)
%          .tuning_param                 - (string) resp table variable name for p-values
%          .tuning_alpha                 - (double) p-value cutoff threshold (default: 0.05)
%          .error_bar_type               - 'sem' (default) or 'std'
%          .region_colors                - N_regions x 3 RGB matrix (optional)
%
% Output:
%   op_out - Output options struct with resolved defaults and figures metadata

%% 1. Validate & Extract Subjects from Input Table
if ~istable(subs)
    error('The second argument (subs) must be a MATLAB table.');
end

if ~ismember('sub', subs.Properties.VariableNames)
    error('The subs table must contain a variable named ''sub''.');
end

% Extract subject IDs as a clean cell array of strings
raw_subs = subs.sub;
if isnumeric(raw_subs)
    sub_list = cellstr(string(raw_subs));
elseif isstring(raw_subs)
    sub_list = cellstr(raw_subs);
elseif iscategorical(raw_subs)
    sub_list = cellstr(string(raw_subs));
elseif iscell(raw_subs)
    sub_list = raw_subs;
else
    sub_list = cellstr(string(raw_subs));
end

% Filter missing/empty entries and keep unique ordering
sub_list = sub_list(~cellfun(@(x) isempty(x) || any(ismissing(x)), sub_list));
sub_list = unique(sub_list, 'stable');
n_subs = length(sub_list);

if n_subs == 0
    error('No valid subject IDs found in subs.sub.');
end

%% 2. Process Brain Regions
if ~ismember('region', resp.Properties.VariableNames)
    resp = define_brain_regions(resp);
end

%% 3. Options & Defaults
if nargin < 3 || isempty(op), op = struct(); end

if ~isfield(op, 'var_to_plot') || isempty(op.var_to_plot)
    op.var_to_plot = 'mean_dist';
end

if ~isfield(op, 'error_bar_type') || isempty(op.error_bar_type)
    op.error_bar_type = 'sem';
end

% Responsivity & Tuning Filter Options
if ~isfield(op, 'analyze_responsive_elcs_only') || isempty(op.analyze_responsive_elcs_only)
    op.analyze_responsive_elcs_only = true;
end

if ~isfield(op, 'analyze_tuned_elcs_only') || isempty(op.analyze_tuned_elcs_only)
    op.analyze_tuned_elcs_only = false;
end

% Validate Tuning Parameter Settings
if op.analyze_tuned_elcs_only
    if ~isfield(op, 'tuning_param') || isempty(op.tuning_param)
        error('op.tuning_param must be specified when op.analyze_tuned_elcs_only is true.');
    end
    if ~ismember(op.tuning_param, resp.Properties.VariableNames)
        error('op.tuning_param ''%s'' was not found in resp table.', op.tuning_param);
    end
    if ~isfield(op, 'tuning_alpha') || isempty(op.tuning_alpha)
        op.tuning_alpha = 0.05;
    end
end

% Determine Responsivity Variable Name in Table
resp_var_name = '';
if op.analyze_responsive_elcs_only
    if ismember('rspv', resp.Properties.VariableNames)
        resp_var_name = 'rspv';
    elseif ismember('rsvp', resp.Properties.VariableNames)
        resp_var_name = 'rsvp';
    else
        warning('Neither ''rspv'' nor ''rsvp'' variable found in resp table. Skipping responsivity filter.');
    end
end

% Find a valid non-empty cnob entry to extract defaults if needed
valid_idx = find(~cellfun(@isempty, resp.cnob), 1, 'first');
if isempty(valid_idx)
    error('No non-empty cnob data found in resp table.');
end

sample_cnob = resp.cnob{valid_idx};

% Resolve Epochs to Plot
if ~isfield(op, 'epochs_to_analyze') || isempty(op.epochs_to_analyze)
    epochs_to_plot = sample_cnob.Properties.RowNames(:)';
else
    epochs_to_plot = op.epochs_to_analyze(:)';
end

% Resolve Conditions to Plot & Order
sample_trialconds = sample_cnob.trialconds{1};
if ~isfield(op, 'sort_cond_vals') || isempty(op.sort_cond_vals)
    conds_to_plot = sample_trialconds.Properties.RowNames(:)';
else
    conds_to_plot = op.sort_cond_vals(:)';
end

% Get Unique Regions (excluding missing/empty)
raw_regions = resp.region;
if iscell(raw_regions)
    valid_mask = ~cellfun(@(x) isempty(x) || any(ismissing(x)), raw_regions);
    unique_regions = unique(raw_regions(valid_mask), 'stable');
else
    unique_regions = unique(raw_regions(~isnan(raw_regions)), 'stable');
end

n_regions = length(unique_regions);
n_conds = length(conds_to_plot);
n_epochs = length(epochs_to_plot);

% Setup Region Colors
if ~isfield(op, 'region_colors') || isempty(op.region_colors)
    colors = lines(n_regions);
else
    colors = op.region_colors;
end

%% 4. Monitor Setup for Alternating Display
monitors = get(0, 'MonitorPositions');
n_monitors = size(monitors, 1);

fig_handles = gobjects(n_epochs, 1);

%% 5. Main Plotting Loop (One Figure per Epoch)
for iepoch = 1:n_epochs
    this_epoch = epochs_to_plot{iepoch};
    
    % Select target monitor for alternating display
    target_mon = mod(iepoch - 1, n_monitors) + 1;
    mon_pos = monitors(target_mon, :);
    
    % Create maximized figure on target monitor
    fig_handles(iepoch) = figure('Name', sprintf('Crossnobis [%s] - Epoch: %s', op.var_to_plot, this_epoch), ...
        'Color', 'w', 'Position', [mon_pos(1)+50, mon_pos(2)+50, mon_pos(3)-100, mon_pos(4)-100]);
    
    try
        fig_handles(iepoch).WindowState = 'maximized';
    catch
        % Fallback if WindowState is unsupported
    end
    
    % Total subplots = N subjects + 1 for combined 'all_subjects'
    n_subplots = n_subs + 1;
    n_cols = ceil(sqrt(n_subplots));
    n_rows = ceil(n_subplots / n_cols);
    
    for isub = 1:n_subplots
        ax = subplot(n_rows, n_cols, isub);
        hold(ax, 'on');
        
        % Filter Subject Data
        if isub <= n_subs
            sub_id = sub_list{isub};
            sub_mask = strcmp(string(resp.sub), string(sub_id));
            plot_title = sprintf('Sub-%s', sub_id);
        else
            sub_mask = ismember(string(resp.sub), string(sub_list)); % Combine selected subjects
            plot_title = 'All Subjects';
        end
        
        % Initialize matrices: [n_conds x n_regions]
        bar_data = nan(n_conds, n_regions);
        err_data = nan(n_conds, n_regions);
        
        for icond = 1:n_conds
            this_cond = conds_to_plot{icond};
            
            for ireg = 1:n_regions
                this_reg = unique_regions{ireg};
                
                % Base electrode mask (subject & region)
                elc_mask = sub_mask & strcmp(resp.region, this_reg);
                
                % Exclude bad electrodes
                if ismember('bad_elc', resp.Properties.VariableNames)
                    elc_mask = elc_mask & ~resp.bad_elc;
                end
                
                % Filter by general task responsivity
                if op.analyze_responsive_elcs_only && ~isempty(resp_var_name)
                    elc_mask = elc_mask & logical(resp.(resp_var_name));
                end
                
                % Filter by feature tuning
                if op.analyze_tuned_elcs_only
                    p_vals = resp.(op.tuning_param);
                    if iscell(p_vals)
                        p_vals = cell2mat(p_vals);
                    end
                    elc_mask = elc_mask & (p_vals < op.tuning_alpha);
                end
                
                elc_indices = find(elc_mask);
                vals = [];
                
                for idx = 1:length(elc_indices)
                    ichan = elc_indices(idx);
                    if isempty(resp.cnob{ichan}), continue; end
                    
                    cnob_tbl = resp.cnob{ichan};
                    if ~ismember(this_epoch, cnob_tbl.Properties.RowNames), continue; end
                    
                    trialconds_tbl = cnob_tbl.trialconds{this_epoch};
                    if ~ismember(this_cond, trialconds_tbl.Properties.RowNames), continue; end
                    
                    v = trialconds_tbl{this_cond, op.var_to_plot};
                    if ~isnan(v)
                        vals(end+1) = v; %#ok<AGROW>
                    end
                end
                
                % Compute Mean and Error
                if ~isempty(vals)
                    bar_data(icond, ireg) = mean(vals);
                    if strcmp(op.error_bar_type, 'sem')
                        err_data(icond, ireg) = std(vals) / sqrt(length(vals));
                    else
                        err_data(icond, ireg) = std(vals);
                    end
                end
            end
        end
        
        % Plot Grouped Bars (X-axis = Conditions, Bars = Regions)
        hb = bar(ax, 1:n_conds, bar_data, 'grouped', 'EdgeColor', 'none');
        
        % Apply Region Colors to Bars & Compute Errorbar X-Positions
        for ireg = 1:n_regions
            hb(ireg).FaceColor = colors(ireg, :);
            
            x_coords = hb(ireg).XEndPoints;
            y_coords = bar_data(:, ireg)';
            e_coords = err_data(:, ireg)';
            
            valid_e = ~isnan(y_coords) & ~isnan(e_coords);
            if any(valid_e)
                errorbar(ax, x_coords(valid_e), y_coords(valid_e), e_coords(valid_e), ...
                    'k.', 'LineWidth', 1.2, 'CapSize', 6);
            end
        end
        
        % Formatting Subplot
        title(ax, plot_title, 'FontSize', 11, 'FontWeight', 'bold');
        set(ax, 'XTick', 1:n_conds, 'XTickLabel', conds_to_plot, 'FontSize', 9, 'TickDir', 'out');
        grid(ax, 'on');
        ax.GridAlpha = 0.15;
        
        if isub == 1 || mod(isub-1, n_cols) == 0
            ylabel(ax, sprintf('%s (%s)', strrep(op.var_to_plot, '_', ' '), upper(op.error_bar_type)));
        end
        
        if isub == n_subplots
            legend(ax, hb, unique_regions, 'Location', 'northeastoutside', 'Interpreter', 'none');
        end
        
        hold(ax, 'off');
    end
    
    % Main Figure Title Annotation
    filter_str = '';
    if op.analyze_responsive_elcs_only
        filter_str = [filter_str, ' | Responsive Only']; 
    end
    if op.analyze_tuned_elcs_only
        filter_str = [filter_str, sprintf(' | Tuned Only (%s < %.3f)', op.tuning_param, op.tuning_alpha)]; 
    end
    
    sgtitle(sprintf('Epoch: %s | Variable: %s%s', this_epoch, op.var_to_plot, filter_str), ...
        'FontSize', 13, 'FontWeight', 'bold');
end

op.fig_handles = fig_handles;
op_out = op;
end