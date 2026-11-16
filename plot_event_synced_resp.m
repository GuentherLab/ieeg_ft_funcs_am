function [region_synced, op_out] = plot_event_synced_resp(resp, subs, op)
% PLOT_EVENT_SYNCED_RESP Plots event-synced electrode responses grouped by brain region.
%
% Syntax:
%   [region_synced, op_out] = plot_event_synced_resp(resp, subs, op)
%
% Inputs:
%   resp - Table with 1 row per electrode containing ephys event responses 
%          (e.g., in resp.disf) and anatomical/responsivity metadata.
%   subs - Subject table containing event timing metadata tables (e.g., in subs.disf).
%   op   - Options struct with configurable parameters:
%          .resptab_varname          - Variable in resp & subs holding event data (default 'events')
%          .trace_spec               - Table/struct/cell specifying traces to plot with fields/cols:
%                                      .name  - Trace label override (if empty, auto-generated from spec)
%                                      .spec  - Cell array of key-value criteria, e.g. {{'type','speech_epoch'},{'delay',{0 228}}}
%                                      .color - Custom RGB/[R G B] or string color (if empty, uses trace_colors)
%          .trace_colors             - Default colormap string or N x 3 RGB matrix (default 'lines')
%          .regions_to_plot          - Cell array of region names to plot (default {} -> all regions)
%          .analyze_responsive_elcs_only - 1 to filter resp.rspv==1, 0 otherwise (default 1)
%          .analyze_tuned_elcs_only      - 1 to filter by tuning_param < tuning_alpha, 0 otherwise (default 0)
%          .tuning_param             - Column name in resp for tuning filter
%          .tuning_alpha             - Alpha threshold for tuning filter (default 0.05)
%          .smooth_windowsize_sec    - Moving average window size in seconds (default 0 = no smoothing)
%          .smooth_method            - Smoothing method for smoothdata (default 'gaussian')
%          .sync_window_sec          - [pre_sec post_sec] window around event (default [] = full saved window)
%          .min_events_per_elc       - Minimum occurrences of an event per elc to include (default 1)
%          .min_elcs_per_region      - Minimum electrodes per region to plot (default 1)
%          .errorbar_alpha           - Transparency for SEM shading (default 0.20)
%          .newfig                   - 1 to create a new figure, 0 to plot in existing (default 1)
%          .trace_width              - Line width for mean traces (default 1.5)
%          .row_col_ratio            - Desired aspect ratio for subplot layout (default 4/3)
%
% Outputs:
%   region_synced - Table with 1 row per region containing nested tables of sorted
%                   response means, SEM, and event/electrode counts per trace.
%   op_out        - Updated options struct.

%% 1. Default Parameters
if nargin < 3 || isempty(op)
    op = struct();
end

field_default('op', 'resptab_varname', 'events');
field_default('op', 'regions_to_plot', {});
field_default('op', 'analyze_responsive_elcs_only', 1);
field_default('op', 'analyze_tuned_elcs_only', 0);
field_default('op', 'tuning_alpha', 0.05);
field_default('op', 'smooth_windowsize_sec', 0);
field_default('op', 'smooth_method', 'gaussian');
field_default('op', 'sync_window_sec', []); % Default to [] (full saved window)
field_default('op', 'min_events_per_elc', 1);
field_default('op', 'min_elcs_per_region', 1);

if isfield(op, 'event_colors') && ~isfield(op, 'trace_colors')
    op.trace_colors = op.event_colors;
end
field_default('op', 'trace_colors', 'lines');

field_default('op', 'errorbar_alpha', 0.20);
field_default('op', 'newfig', 1);
field_default('op', 'trace_width', 1.5);
field_default('op', 'row_col_ratio', 4/3);
field_default('op', 'resp_signal', 'hg');

varname = op.resptab_varname;

%% 2. Validate Response & Event Table Variables
if ~ismember(varname, resp.Properties.VariableNames)
    error('plot_event_synced_resp:MissingRespVar', ...
        'Variable "%s" not found in resp table.', varname);
end

if ~ismember(varname, subs.Properties.VariableNames)
    error('plot_event_synced_resp:MissingSubsVar', ...
        'Variable "%s" not found in subs table.', varname);
end

%% 3. Filter Electrodes (Responsivity & Tuning)
if op.analyze_responsive_elcs_only
    if ismember('rspv', resp.Properties.VariableNames)
        resp = resp(logical(resp.rspv), :);
    else
        warning('plot_event_synced_resp:NoRspvColumn', ...
            'op.analyze_responsive_elcs_only==1 but ''rspv'' column missing in resp table. Skipping.');
    end
end

if op.analyze_tuned_elcs_only
    if ~isfield(op, 'tuning_param') || isempty(op.tuning_param)
        error('plot_event_synced_resp:MissingTuningParam', ...
            'op.analyze_tuned_elcs_only==1 but op.tuning_param was not specified.');
    end
    
    if iscell(op.tuning_param)
        t_col = op.tuning_param{1};
        t_idx = op.tuning_param{2};
        tuned_mask = resp{:, t_col}(:, t_idx) < op.tuning_alpha;
    else
        t_col = op.tuning_param;
        tuned_mask = resp{:, t_col} < op.tuning_alpha;
    end
    resp = resp(tuned_mask, :);
end

nelcs = height(resp);
if nelcs == 0
    error('plot_event_synced_resp:NoElectrodes', 'No electrodes passed filtering criteria.');
end

%% 4. Define Brain Regions
[resp, op] = define_brain_regions(resp, op);

%% 5. Parse and Standardize Trace Specifications (op.trace_spec)
if ~isfield(op, 'trace_spec') || isempty(op.trace_spec)
    all_avail_types = {};
    for isub = 1:height(subs)
        ev_tbl = subs.(varname){isub};
        if ~isempty(ev_tbl) && istable(ev_tbl) && ismember('type', ev_tbl.Properties.VariableNames)
            all_avail_types = [all_avail_types; cellstr(string(ev_tbl.type))]; %#ok<AGROW>
        end
    end
    all_avail_types = unique(all_avail_types, 'stable');
    
    n_types = length(all_avail_types);
    names_col = cell(n_types, 1);
    specs_col = cell(n_types, 1);
    colors_col = cell(n_types, 1);
    for k = 1:n_types
        specs_col{k} = {{'type', all_avail_types{k}}};
    end
    op.trace_spec = table(names_col, specs_col, colors_col, 'VariableNames', {'name', 'spec', 'color'});
end

if istable(op.trace_spec)
    raw_names  = op.trace_spec.name;
    raw_specs  = op.trace_spec.spec;
    raw_colors = op.trace_spec.color;
elseif isstruct(op.trace_spec)
    raw_names  = {op.trace_spec.name}';
    raw_specs  = {op.trace_spec.spec}';
    raw_colors = {op.trace_spec.color}';
elseif iscell(op.trace_spec)
    raw_names  = op.trace_spec(:, 1);
    raw_specs  = op.trace_spec(:, 2);
    raw_colors = op.trace_spec(:, 3);
else
    error('plot_event_synced_resp:InvalidTraceSpec', 'op.trace_spec must be a table, struct, or cell array.');
end

ntraces = length(raw_specs);
trace_base_names = cell(ntraces, 1);

for itrace = 1:ntraces
    user_name = raw_names{itrace};
    if iscell(user_name) && ~isempty(user_name)
        user_name = user_name{1};
    end
    
    if ~isempty(user_name) && (ischar(user_name) || isstring(user_name)) && strlength(string(user_name)) > 0
        trace_base_names{itrace} = char(user_name);
    else
        % Auto-generate trace name from spec (joining multiple values with '/')
        this_spec = raw_specs{itrace};
        name_parts = {};
        for icond = 1:length(this_spec)
            cond_pair = this_spec{icond};
            col_k = cond_pair{1};
            val_k = cond_pair{2};
            
            val_str = strjoin(string(val_k), '/');
            if strcmpi(col_k, 'type')
                name_parts{end+1} = char(val_str); %#ok<AGROW>
            else
                name_parts{end+1} = sprintf('%s=%s', col_k, val_str); %#ok<AGROW>
            end
        end
        trace_base_names{itrace} = strjoin(name_parts, ', ');
    end
end

%% 6. Determine Saved Window and Calculate Relative Time Vector
sample_mat = [];
for ielc = 1:nelcs
    if ~isempty(resp.(varname){ielc})
        sample_mat = resp.(varname){ielc};
        break;
    end
end

if isempty(sample_mat)
    error('plot_event_synced_resp:EmptyData', 'No valid event response matrices found in resp.%s', varname);
end

ntimepoints_full = size(sample_mat, 2);

if isfield(op, 'sync_window_sec_saved') && ~isempty(op.sync_window_sec_saved)
    full_sync_win = op.sync_window_sec_saved;
else
    full_sync_win = [3 3];
end

times_full = linspace(-full_sync_win(1), full_sync_win(2), ntimepoints_full);
dt = (full_sync_win(1) + full_sync_win(2)) / (ntimepoints_full - 1);

if ~isempty(op.sync_window_sec)
    pre_req  = op.sync_window_sec(1);
    post_req = op.sync_window_sec(2);
    col_mask = (times_full >= -pre_req - 1e-5) & (times_full <= post_req + 1e-5);
    times_rel = times_full(col_mask);
else
    col_mask = true(1, ntimepoints_full);
    times_rel = times_full;
end

%% 7. Filter Regions to Plot
if ~isempty(op.regions_to_plot)
    if ~iscell(op.regions_to_plot)
        op.regions_to_plot = cellstr(string(op.regions_to_plot));
    end
    
    region_names_all = op.regiondef.region;
    valid_indices = [];
    for i = 1:length(op.regions_to_plot)
        req_reg = op.regions_to_plot{i};
        idx = find(strcmp(region_names_all, req_reg), 1);
        if isempty(idx) || isempty(op.regiondef.n_elcs(idx)) || op.regiondef.n_elcs(idx) == 0
            fprintf(1, '    Warning: Region "%s" has no qualifying electrodes.\n', req_reg);
        else
            valid_indices(end+1) = idx; %#ok<AGROW>
        end
    end
    
    if isempty(valid_indices)
        error('plot_event_synced_resp:NoValidRegions', 'None of the specified regions contain valid electrodes.');
    end
    
    op.regiondef = op.regiondef(valid_indices, :);
    op.nregions  = height(op.regiondef);
end

nregions = op.nregions;

%% 8. Configure Tiled Figure Layout
if op.newfig
    figure('Color', 'w', 'WindowState', 'maximized');
end

ntiles = nregions;
r = 1:ntiles;
c = ceil(ntiles ./ r);
[~, idx] = min(abs(c ./ r - 1/op.row_col_ratio));
op.n_plot_rows = r(idx);
op.n_plot_cols = c(idx);

tiledlayout(op.n_plot_rows, op.n_plot_cols, 'Padding', 'none', 'TileSpacing', 'tight');

% Construct Figure Title
nsubs_total = length(unique(resp.sub));
if op.analyze_responsive_elcs_only
    resp_str = 'Responsive Elcs Only';
else
    resp_str = 'All Elcs';
end

if op.analyze_tuned_elcs_only
    if iscell(op.tuning_param)
        t_param_str = op.tuning_param{1};
    else
        t_param_str = op.tuning_param;
    end
    tune_str = sprintf(', Tuned to %s (p<%g)', t_param_str, op.tuning_alpha);
else
    tune_str = '';
end

fig_title = sprintf('Event-Synced %s Power (%s%s) | Subs: %d, Elcs: %d', ...
    upper(op.resp_signal), resp_str, tune_str, nsubs_total, nelcs);
sgtitle(fig_title, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

%% 9. Color Palette Resolution
if ischar(op.trace_colors) || isstring(op.trace_colors)
    default_cmap = feval(char(op.trace_colors), ntraces);
else
    default_cmap = op.trace_colors;
end

trace_colors_final = zeros(ntraces, 3);
for itrace = 1:ntraces
    this_c = raw_colors{itrace};
    if iscell(this_c) && ~isempty(this_c)
        this_c = this_c{1};
    end
    
    if ~isempty(this_c)
        if isnumeric(this_c) && numel(this_c) == 3
            trace_colors_final(itrace, :) = this_c(:)';
        elseif ischar(this_c) || isstring(this_c)
            switch lower(char(this_c))
                case {'r','red'},     trace_colors_final(itrace,:) = [1 0 0];
                case {'g','green'},   trace_colors_final(itrace,:) = [0 0.8 0];
                case {'b','blue'},    trace_colors_final(itrace,:) = [0 0 1];
                case {'k','black'},   trace_colors_final(itrace,:) = [0 0 0];
                case {'m','magenta'}, trace_colors_final(itrace,:) = [1 0 1];
                case {'c','cyan'},    trace_colors_final(itrace,:) = [0 1 1];
                case {'y','yellow'},  trace_colors_final(itrace,:) = [0.9 0.8 0];
                otherwise,            trace_colors_final(itrace,:) = default_cmap(mod(itrace-1, size(default_cmap,1))+1, :);
            end
        else
            trace_colors_final(itrace, :) = default_cmap(mod(itrace-1, size(default_cmap,1))+1, :);
        end
    else
        trace_colors_final(itrace, :) = default_cmap(mod(itrace-1, size(default_cmap,1))+1, :);
    end
end

%% 10. Process & Plot Each Region
region_synced_rows = cell(nregions, 1);
region_names = op.regiondef.region;
nsubs_total_subs = height(subs);

for iregion = 1:nregions
    this_region = region_names{iregion};
    n_elcs_this_region = op.regiondef.n_elcs(iregion);
    hax = nexttile;
    
    elc_indices = find(strcmp(resp.region, this_region));
    nelcs_reg = length(elc_indices);
    
    trace_labels_col = trace_base_names;
    mean_cell = cell(ntraces, 1);
    sem_cell  = cell(ntraces, 1);
    nelc_vec  = zeros(ntraces, 1);
    nevt_vec  = zeros(ntraces, 1);
    
    if nelcs_reg >= op.min_elcs_per_region
        hold(hax, 'on');
        
        reg_ymax = -Inf;
        reg_ymin = Inf;
        
        for itrace = 1:ntraces
            this_spec = raw_specs{itrace};
            base_label = trace_base_names{itrace};
            this_color = trace_colors_final(itrace, :);
            
            elc_means = [];
            
            % Track unique matching event rows per subject for this region & trace
            sub_evt_included = cell(nsubs_total_subs, 1);
            for isub = 1:nsubs_total_subs
                if ~isempty(subs.(varname){isub}) && istable(subs.(varname){isub})
                    sub_evt_included{isub} = false(height(subs.(varname){isub}), 1);
                end
            end
            
            for i = 1:nelcs_reg
                ielc = elc_indices(i);
                sub_name = resp.sub{ielc};
                subind = find(string(subs.sub) == string(sub_name), 1);
                
                if isempty(subind), continue; end
                
                events_tbl = subs.(varname){subind};
                elc_mat    = resp.(varname){ielc};
                
                if isempty(events_tbl) || isempty(elc_mat), continue; end
                
                % Multi-column & multi-value filter mask
                evt_mask = true(height(events_tbl), 1);
                for icond = 1:length(this_spec)
                    col_k = this_spec{icond}{1};
                    val_k = this_spec{icond}{2};
                    
                    if ismember(col_k, events_tbl.Properties.VariableNames)
                        col_vals = events_tbl.(col_k);
                        if isnumeric(col_vals)
                            if iscell(val_k)
                                val_vec = cell2mat(val_k);
                            else
                                val_vec = val_k;
                            end
                            evt_mask = evt_mask & ismember(col_vals, val_vec);
                        else
                            evt_mask = evt_mask & ismember(string(col_vals), string(val_k));
                        end
                    else
                        evt_mask = false(height(events_tbl), 1);
                        break;
                    end
                end
                
                % Filter to rows with valid (non-NaN) data in elc_mat for this electrode
                valid_data_mask = any(~isnan(elc_mat(:, col_mask)), 2);
                elc_evt_mask = evt_mask & valid_data_mask;
                
                n_evts_this = nnz(elc_evt_mask);
                if n_evts_this >= op.min_events_per_elc
                    mat_sub = elc_mat(elc_evt_mask, col_mask);
                    elc_means = [elc_means; nanmean(mat_sub, 1)]; %#ok<AGROW>
                    
                    % Record event rows covered by at least 1 valid electrode
                    if ~isempty(sub_evt_included{subind})
                        sub_evt_included{subind} = sub_evt_included{subind} | elc_evt_mask;
                    end
                end
            end
            
            % Count total unique events meeting spec across subjects with >=1 valid electrode response
            total_events_type = 0;
            for isub = 1:nsubs_total_subs
                if ~isempty(sub_evt_included{isub})
                    total_events_type = total_events_type + nnz(sub_evt_included{isub});
                end
            end
            
            disp_label = sprintf('%s (n=%d)', base_label, total_events_type);
            
            if ~isempty(elc_means)
                mean_tc = nanmean(elc_means, 1);
                sem_tc  = nanstd(elc_means, [], 1) ./ sqrt(sum(~isnan(elc_means), 1));
                
                if op.smooth_windowsize_sec > 0
                    smooth_pts = round(op.smooth_windowsize_sec / dt);
                    if smooth_pts > 1
                        mean_tc = smoothdata(mean_tc, 2, op.smooth_method, smooth_pts);
                        sem_tc  = smoothdata(sem_tc,  2, op.smooth_method, smooth_pts);
                    end
                end
                
                mean_cell{itrace} = mean_tc;
                sem_cell{itrace}  = sem_tc;
                nelc_vec(itrace)  = size(elc_means, 1);
                nevt_vec(itrace)  = total_events_type;
                
                low_b  = mean_tc - sem_tc;
                high_b = mean_tc + sem_tc;
                reg_ymax = max([reg_ymax, max(high_b)]);
                reg_ymin = min([reg_ymin, min(low_b)]);
                
                fill_x = [times_rel, fliplr(times_rel)];
                fill_y = [low_b, fliplr(high_b)];
                
                h_fill = fill(hax, fill_x, fill_y, this_color, ...
                    'FaceAlpha', op.errorbar_alpha, 'EdgeColor', 'none');
                set(h_fill, 'HandleVisibility', 'off');
                
                plot(hax, times_rel, mean_tc, 'LineWidth', op.trace_width, ...
                    'Color', this_color, 'DisplayName', disp_label);
            else
                plot(hax, NaN, NaN, 'LineWidth', op.trace_width, ...
                    'Color', this_color, 'DisplayName', disp_label);
            end
        end
        
        xline(hax, 0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'HandleVisibility', 'off');
        yline(hax, 0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
        xlim(hax, [min(times_rel), max(times_rel)]);
        xlabel(hax, 'Time (sec)');
        ylabel(hax, sprintf('normed %s power', op.resp_signal));
        title(hax, sprintf('%s (n=%d elc)', this_region, n_elcs_this_region), ...
            'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
        box(hax, 'off');
        
        if iregion == 1 && isfinite(reg_ymax) && isfinite(reg_ymin)
            y_range = reg_ymax - reg_ymin;
            if y_range == 0, y_range = 1; end
            
            ylim(hax, [reg_ymin - 0.05 * y_range, reg_ymax + 0.35 * y_range]);
            
            h_leg = legend(hax, 'Location', 'northwest', 'Interpreter', 'none', ...
                'FontSize', 8, 'Box', 'on');
            title(h_leg, 'Sync Event Traces', 'FontSize', 8.5, 'FontWeight', 'bold');
        end
        
        hold(hax, 'off');
    else
        axis(hax, 'off');
    end
    
    synced_summary = table(trace_labels_col, mean_cell, sem_cell, nelc_vec, nevt_vec, ...
        'VariableNames', {'trace_label', 'mean', 'sem', 'n_elcs', 'n_events_total'});
    
    region_synced_rows{iregion} = table(categorical({this_region}), n_elcs_this_region, {synced_summary}, ...
        'VariableNames', {'region', 'n_elcs', 'synced_summary'});
end

%% 11. Finalize Output Table
region_synced = vertcat(region_synced_rows{:});
op_out = op;

end

%% Helper Subfunction
function field_default(s_name, f, v)
    if ischar(s_name) || isstring(s_name)
        s_val = evalin('caller', s_name);
        if ~isfield(s_val, f)
            s_val.(f) = v;
            assignin('caller', s_name, s_val);
        end
    end
end