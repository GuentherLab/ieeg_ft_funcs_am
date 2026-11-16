function [resp_out, events_out] = event_synced_responses(resp_in, D, events_in, cfg)
% EVENT_SYNCED_RESPONSES Outputs electrode ephys responses time-locked to events,
% with optional trial-matched baselining.
%
% Syntax:
%   [resp_out, events_out] = event_synced_responses(resp_in, D, events_in, cfg)
%
% Inputs:
%   resp_in    - Table with 1 row per electrode (must contain 'chan' column, 
%                and 'base' column if cfg.do_baselining is true).
%   D          - FieldTrip raw data structure containing ephys data:
%                .trial{1} - channel x time matrix
%                .time{1}  - 1D vector of global time coordinates (GTC, in seconds)
%                .label    - cell array of channel labels
%                .fsample  - sampling rate in Hz (optional)
%   events_in  - Table containing 1 row per event (must contain 'type' and timestamp column).
%   cfg        - Configuration struct with fields:
%                .event_types               - cell array of strings (matches events_in.type)
%                .sync_window_sec           - [pre_sec, post_sec] window around event (default [1 1])
%                .resptab_out_varname       - name of variable added to resp_out (default 'events')
%                .event_time_col            - column name in events_in for timestamp (default 'starts')
%                .do_baselining             - logical (default true); if false, outputs raw responses
%                .trials                    - behavioral trials table (required if do_baselining=true)
%                .epochs                    - epochs table containing 'base' row (required if do_baselining=true)
%                .baseline_method           - 'subtract_then_divide' (default), 'subtract', etc.
%                .max_timecourse_base_ratio - ratio threshold to exclude artifactual windows (default 50)
%
% Outputs:
%   resp_out   - Input table with a new cell column containing a (nevents x ntimepoints)
%                matrix per electrode.
%   events_out - Filtered events table matching cfg.event_types.

%% 1. Configuration & Default Options
if nargin < 4 || isempty(cfg)
    cfg = struct();
end

% Baselining options & defaults
if ~isfield(cfg, 'do_baselining') || isempty(cfg.do_baselining)
    cfg.do_baselining = true;
end

if ~isfield(cfg, 'sync_window_sec') || isempty(cfg.sync_window_sec)
    cfg.sync_window_sec = [1 1];
end

if ~isfield(cfg, 'resptab_out_varname') || isempty(cfg.resptab_out_varname)
    cfg.resptab_out_varname = 'events';
end

if ~isfield(cfg, 'event_types')
    cfg.event_types = {};
end

% Validate baselining inputs when baselining is enabled
if cfg.do_baselining
    if ~isfield(cfg, 'trials') || isempty(cfg.trials) || ~istable(cfg.trials)
        error('event_synced_responses:InvalidTrialsTable', ...
            'cfg.trials must be a valid, non-empty table when cfg.do_baselining is true.');
    end
    if ~isfield(cfg, 'epochs') || isempty(cfg.epochs) || ~istable(cfg.epochs)
        error('event_synced_responses:InvalidEpochsTable', ...
            'cfg.epochs must be a valid, non-empty table when cfg.do_baselining is true.');
    end
    if ~ismember('base', resp_in.Properties.VariableNames)
        error('event_synced_responses:MissingRespBaseColumn', ...
            'resp_in must contain a ''base'' column when cfg.do_baselining is true.');
    end
    
    % Locate 'base' spec in cfg.epochs
    if ismember('base', cfg.epochs.Properties.RowNames)
        base_spec_offset = cfg.epochs.offset{'base'};
    elseif ismember('epoch', cfg.epochs.Properties.VariableNames) && ismember('base', cfg.epochs.epoch)
        base_idx = find(strcmp(cfg.epochs.epoch, 'base'), 1);
        base_spec_offset = cfg.epochs.offset{base_idx};
    else
        error('event_synced_responses:MissingBaseEpoch', ...
            'cfg.epochs table must contain a row or entry named ''base''.');
    end

    if ~isfield(cfg, 'baseline_method') || isempty(cfg.baseline_method)
        cfg.baseline_method = 'subtract_then_divide';
    end
    if ~isfield(cfg, 'max_timecourse_base_ratio') || isempty(cfg.max_timecourse_base_ratio)
        cfg.max_timecourse_base_ratio = 50;
    end
end

%% 2. Extract Sampling Rate and Calculate Timepoint Window Size
if isfield(D, 'fsample') && ~isempty(D.fsample)
    fs = D.fsample;
elseif isfield(D, 'time')
    if iscell(D.time)
        fs = 1 / mean(diff(D.time{1}));
    else
        fs = 1 / mean(diff(D.time));
    end
else
    error('event_synced_responses:MissingFsample', 'Sampling rate (fsample) could not be determined from D.');
end

n_pre  = round(cfg.sync_window_sec(1) * fs);
n_post = round(cfg.sync_window_sec(2) * fs);
ntimepoints = n_pre + n_post + 1;

%% 3. Filter Events Table (events_out)
if ~ismember('type', events_in.Properties.VariableNames)
    error('event_synced_responses:MissingTypeColumn', 'events_in table must contain a ''type'' column.');
end

% Identify event timestamp column in events_in
if isfield(cfg, 'event_time_col') && ismember(cfg.event_time_col, events_in.Properties.VariableNames)
    event_time_col = cfg.event_time_col;
else
    candidate_cols = {'starts', 'onset', 'time', 't_onset', 't_start', 'start'};
    matched_idx = find(ismember(candidate_cols, events_in.Properties.VariableNames), 1);
    if isempty(matched_idx)
        error('event_synced_responses:MissingTimestampColumn', 'Could not find event timestamp column in events_in.');
    end
    event_time_col = candidate_cols{matched_idx};
end

events_type_str = string(events_in.type);
if ~isempty(cfg.event_types)
    target_types_str = string(cfg.event_types);
    valid_event_mask = ismember(events_type_str, target_types_str);
    events_out = events_in(valid_event_mask, :);
else
    events_out = events_in;
end

nevents = height(events_out);

%% 4. Unbox FieldTrip Data Matrix and Global Time Vector
if iscell(D.trial)
    data_mat = D.trial{1};
else
    data_mat = D.trial;
end

if iscell(D.time)
    time_vec = D.time{1};
else
    time_vec = D.time;
end

nelcs = height(resp_in);
d_labels = string(D.label);
resp_col = cell(nelcs, 1);

%% 5. Match Each Event to the Most Recent Preceding Trial Baseline
matched_trial_per_event = nan(nevents, 1);
if cfg.do_baselining && nevents > 0
    ntrials_beh = height(cfg.trials);
    trial_base_offsets = nan(ntrials_beh, 1);
    
    for itr = 1:ntrials_beh
        trial_base_offsets(itr) = parse_epoch_time(base_spec_offset, cfg.trials, itr);
    end
    
    event_times = double(events_out.(event_time_col));
    for ievt = 1:nevents
        t_event = event_times(ievt);
        if isnan(t_event)
            continue;
        end
        
        % The most recent baseline must end before the time-zero of this event
        valid_trials = find(trial_base_offsets < t_event);
        if ~isempty(valid_trials)
            [~, max_rel_idx] = max(trial_base_offsets(valid_trials));
            matched_trial_per_event(ievt) = valid_trials(max_rel_idx);
        end
    end
else
    if nevents > 0
        event_times = double(events_out.(event_time_col));
    end
end

%% 6. Extract and Normalize Event-Synced Response Matrices (nevents x ntimepoints)
if nevents == 0
    warning('event_synced_responses:NoMatchingEvents', ...
        'No events matched specified cfg.event_types. Returning empty matrices of size (0 x %d).', ntimepoints);
    for ielc = 1:nelcs
        resp_col{ielc} = nan(0, ntimepoints);
    end
else
    t_start_global = time_vec(1);
    total_samples  = length(time_vec);

    for ielc = 1:nelcs
        if iscell(resp_in.chan)
            chan_name = string(resp_in.chan{ielc});
        else
            chan_name = string(resp_in.chan(ielc));
        end

        d_chan_idx = find(d_labels == chan_name, 1);
        if isempty(d_chan_idx)
            warning('event_synced_responses:ChannelNotFound', ...
                'Channel "%s" not found in D.label. Filling with NaNs.', chan_name);
            resp_col{ielc} = nan(nevents, ntimepoints);
            continue;
        end

        elc_mat = nan(nevents, ntimepoints);
        raw_chan_signal = data_mat(d_chan_idx, :);

        for ievt = 1:nevents
            t0 = event_times(ievt);
            if isnan(t0)
                continue;
            end

            % Find center sample index corresponding to event onset
            idx_center = round((t0 - t_start_global) * fs) + 1;
            idx_start  = idx_center - n_pre;
            idx_end    = idx_start + ntimepoints - 1;

            % Extract raw signal slice with boundary checking
            v_start = max(1, idx_start);
            v_end   = min(total_samples, idx_end);

            raw_win = nan(1, ntimepoints);
            if v_start <= v_end
                out_start = (v_start - idx_start) + 1;
                out_end   = out_start + (v_end - v_start);
                raw_win(out_start:out_end) = raw_chan_signal(v_start:v_end);
            end

            % Perform baseline normalization if enabled
            if cfg.do_baselining
                matched_trial_idx = matched_trial_per_event(ievt);
                if isnan(matched_trial_idx)
                    elc_mat(ievt, :) = nan(1, ntimepoints);
                else
                    base_val = resp_in.base{ielc}(matched_trial_idx);
                    if isnan(base_val) || base_val == 0
                        elc_mat(ievt, :) = nan(1, ntimepoints);
                    else
                        normed_win = do_baselining_calc(raw_win, base_val, cfg.baseline_method);
                        if ~isempty(cfg.max_timecourse_base_ratio) && isfinite(cfg.max_timecourse_base_ratio) && ...
                                max(normed_win) > cfg.max_timecourse_base_ratio
                            elc_mat(ievt, :) = nan(1, ntimepoints);
                        else
                            elc_mat(ievt, :) = normed_win;
                        end
                    end
                end
            else
                elc_mat(ievt, :) = raw_win;
            end
        end

        resp_col{ielc} = elc_mat;
    end
end

%% 7. Assign Column to Output Table
resp_out = resp_in;
out_varname = cfg.resptab_out_varname;
resp_out.(out_varname) = resp_col;

end

%% Helper Functions

function t = parse_epoch_time(spec, trial_table, itrial)
    if iscell(spec)
        if length(spec) == 1 && iscell(spec{1})
            spec = spec{1};
        end
        evt_name = spec{1};
        if length(spec) >= 2
            shift = spec{2};
        else
            shift = 0;
        end
    else
        evt_name = spec;
        shift = 0;
    end
    t = trial_table{itrial, evt_name} + shift;
end

function normed_response = do_baselining_calc(response, baseval, method)
    switch method
        case 'subtract'
            normed_response = response - baseval;
        case 'subtract_then_divide'
            normed_response = (response - baseval) / baseval;
        case 'divide_then_subtract'
            normed_response = (response / baseval) - 1;
        case 'divide'
            normed_response = response / baseval;
        otherwise
            error('event_synced_responses:UnknownBaselineMethod', ...
                'Unrecognized baseline method: %s', method);
    end
end