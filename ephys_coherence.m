function coh = ephys_coherence(subs, op)
% EPHYS_COHERENCE - Computes pairwise inter-regional trial-aligned coherence across ephys channels.
%
% Usage:
%   coh = ephys_coherence(subs, op)
%
% Inputs:
%   subs - Table containing subject path definitions in column 'paths' (nsubjects x 1 cell array).
%          Each cell contains a struct with subject-specific paths/variables:
%            .fieldtrip  : Path to .mat file containing FieldTrip continuous data struct
%            .trials     : (Optional) Path to .tsv table OR pre-loaded table containing trial timestamps
%            .artifact   : Path or cell array of paths to manual artifact .tsv files
%            .electrodes : Path to .tsv table listing electrode metadata
%            .resp       : Path to .mat file, .tsv table, OR pre-loaded table containing 'resp' metadata
%          If paths.trials is missing/empty, the function uses subs.trials{isub}.
%
%   op   - Options structure with fields:
%            .sync_events        : Table with rows/rownames as event names and columns 'start', 'end' (sec)
%            .freqs              : Numerical vector of target frequencies (e.g., [20, 100])
%            .regions_to_analyze : Cell array of region strings to pair inter-regionally (e.g., {'IFG/IFS','STN','Thal'})
%            .sort_cond          : String specifying the column in trials table to group trials by
%            .sort_cond_vals     : (Optional) Cell array/vector specifying exact condition values and order
%            .resp_vars_to_copy  : Cell array of column names to copy from resp table into output coh table
%            .savefile           : Filepath string for saving output MAT-file
%            .coh_measure        : Coherence metric - 'mag_sq' (default) or 'imag'
%            .buffer_window_sec  : Buffer duration in seconds for wavelet edge artifact suppression (default 1.0)
%            .fieldtrip_var_name : Name of variable in fieldtrip file (default 'D_ref')
%            .show_progress      : (Optional) Logical flag (true/false) to display command line progress per subject
%
% Outputs:
%   coh  - Table containing 1 row per inter-regional channel pair across all subjects with nested coherence data.

%% =======================================================================
%% 1. DEFAULT OPTION ASSIGNMENTS & INITIAL VALIDATION
%% =======================================================================
if nargin < 2 || isempty(op)
    op = struct();
end

% Set default parameter values if omitted by user
if ~isfield(op, 'fieldtrip_var_name') || isempty(op.fieldtrip_var_name)
    op.fieldtrip_var_name = 'D_ref';
end

if ~isfield(op, 'coh_measure') || isempty(op.coh_measure)
    % Default to magnitude squared coherence as used in Kingyon et al. 2015
    % Note: Kingyon et al. (2015) used magnitude-squared coherence calculated via short-time FFT.
    op.coh_measure = 'mag_sq';
end

if ~isfield(op, 'buffer_window_sec') || isempty(op.buffer_window_sec)
    % 1-second edge buffer added to start/end before wavelet convolution to eliminate edge artifacts
    op.buffer_window_sec = 1.0;
end

if ~isfield(op, 'resp_vars_to_copy')
    op.resp_vars_to_copy = {};
end

if ~isfield(op, 'sort_cond_vals')
    op.sort_cond_vals = {};
end

if ~isfield(op, 'show_progress') || isempty(op.show_progress)
    op.show_progress = false;
end

% Determine number of subjects
nsubjects = height(subs);
master_pair_tables = cell(nsubjects, 1);

%% =======================================================================
%% 2. SUBJECT-LEVEL PROCESSING LOOP
%% =======================================================================
for isub = 1:nsubjects
    % Determine subject identifier label
    if ismember('sub', subs.Properties.VariableNames)
        sub_name = subs.sub{isub};
    elseif ismember('name', subs.Properties.VariableNames)
        sub_name = subs.name{isub};
    else
        sub_name = sprintf('sub-%02d', isub);
    end
    
    % Display progress if requested
    if op.show_progress
        fprintf('Processing subject %d/%d: %s...\n', isub, nsubjects, sub_name);
    end
    
    paths = subs.paths{isub};
    
    %% --- 2A. Load FieldTrip Continuous Data Structure ---
    if ~exist(paths.fieldtrip, 'file')
        error('ephys_coherence:FileNotFound', 'Fieldtrip file not found: %s', paths.fieldtrip);
    end
    
    ft_file_info = whos('-file', paths.fieldtrip);
    ft_var_names = {ft_file_info.name};
    
    if ~ismember(op.fieldtrip_var_name, ft_var_names)
        error('ephys_coherence:VarNotFound', ...
            'Expected variable ''%s'' in file ''%s'', but found variables: %s', ...
            op.fieldtrip_var_name, paths.fieldtrip, strjoin(ft_var_names, ', '));
    end
    
    loaded_ft = load(paths.fieldtrip, op.fieldtrip_var_name);
    D = loaded_ft.(op.fieldtrip_var_name);
    
    % Extract sampling frequency and time vector
    time_vec = D.time{1};
    if isfield(D, 'fsample') && ~isempty(D.fsample)
        fs = D.fsample;
    else
        fs = 1 / mean(diff(time_vec));
    end
    
    %% --- 2B. Load & Process Electrodes Table ---
    if istable(paths.electrodes)
        electrodes = paths.electrodes;
    else
        electrodes = readtable(paths.electrodes, 'FileType', 'text', 'Delimiter', '\t');
    end
    
    % Define anatomical brain regions using define_brain_regions
    [electrodes, ~] = define_brain_regions(electrodes);
    
    % Map FieldTrip channel labels to electrodes.region
    nchans = length(D.label);
    chan_regions = cell(nchans, 1);
    
    for ichan = 1:nchans
        chan_name = D.label{ichan};
        
        % Identify matching electrode row via electrode_label or name
        if isfield(D, 'electrode_label') && length(D.electrode_label) >= ichan && ~isempty(D.electrode_label{ichan})
            match_idx = find(strcmp(string(electrodes.name), string(D.electrode_label{ichan})), 1);
        else
            match_idx = find(strcmp(string(electrodes.name), string(chan_name)), 1);
        end
        
        if ~isempty(match_idx) && ismember('region', electrodes.Properties.VariableNames) && ...
                ~isempty(electrodes.region{match_idx}) && ~all(ismissing(string(electrodes.region{match_idx})))
            chan_regions{ichan} = char(electrodes.region{match_idx});
        else
            chan_regions{ichan} = 'unknown';
        end
    end
    
    %% --- 2C. Load & Stack Artifact Tables ---
    artifact_master = table();
    if isfield(paths, 'artifact') && ~isempty(paths.artifact)
        art_files = paths.artifact;
        if ischar(art_files) || isstring(art_files)
            art_files = cellstr(art_files);
        end
        
        for iart = 1:length(art_files)
            if exist(art_files{iart}, 'file')
                curr_art = readtable(art_files{iart}, 'FileType', 'text', 'Delimiter', '\t');
                artifact_master = [artifact_master; curr_art]; %#ok<AGROW>
            end
        end
    end
    
    %% --- 2D. Load Trials Table ---
    if isfield(paths, 'trials') && ~isempty(paths.trials)
        if istable(paths.trials)
            trials_tbl = paths.trials;
        else
            trials_tbl = readtable(paths.trials, 'FileType', 'text', 'Delimiter', '\t');
        end
    elseif ismember('trials', subs.Properties.VariableNames) && ~isempty(subs.trials{isub})
        if istable(subs.trials{isub})
            trials_tbl = subs.trials{isub};
        elseif iscell(subs.trials{isub})
            trials_tbl = subs.trials{isub}{1};
        else
            error('ephys_coherence:InvalidTrials', 'subs.trials{%d} must contain a valid MATLAB table.', isub);
        end
    else
        error('ephys_coherence:TrialsNotFound', 'No trials table found in paths.trials or subs.trials for subject %d.', isub);
    end
    
    %% --- 2E. Load Channel Response ('resp') Table ---
    resp_tbl = table();
    if isfield(paths, 'resp') && ~isempty(paths.resp)
        if istable(paths.resp)
            resp_tbl = paths.resp;
        elseif ischar(paths.resp) || isstring(paths.resp)
            resp_path = char(paths.resp);
            if exist(resp_path, 'file')
                [~, ~, ext] = fileparts(resp_path);
                if strcmpi(ext, '.mat')
                    resp_data = load(resp_path);
                    if isfield(resp_data, 'resp')
                        resp_tbl = resp_data.resp;
                    else
                        fnames = fieldnames(resp_data);
                        if ~isempty(fnames)
                            resp_tbl = resp_data.(fnames{1});
                        end
                    end
                elseif strcmpi(ext, '.tsv') || strcmpi(ext, '.txt') || strcmpi(ext, '.csv')
                    resp_tbl = readtable(resp_path, 'FileType', 'text');
                end
            end
        end
    elseif ismember('resp', subs.Properties.VariableNames) && ~isempty(subs.resp{isub})
        if istable(subs.resp{isub})
            resp_tbl = subs.resp{isub};
        elseif ischar(subs.resp{isub}) || isstring(subs.resp{isub})
            resp_path = char(subs.resp{isub});
            if exist(resp_path, 'file')
                resp_data = load(resp_path);
                if isfield(resp_data, 'resp')
                    resp_tbl = resp_data.resp;
                end
            end
        end
    end
    
    %% --- 2F. Select Channels & Construct Inter-Regional Pairs ---
    target_regions = op.regions_to_analyze;
    valid_chan_mask = ismember(string(chan_regions), string(target_regions));
    valid_chan_indices = find(valid_chan_mask);
    
    % Build pairwise inter-regional combinations (exclude intra-regional pairs)
    pair_c1 = [];
    pair_c2 = [];
    
    for i1 = 1:length(valid_chan_indices)
        idx1 = valid_chan_indices(i1);
        reg1 = chan_regions{idx1};
        
        for i2 = (i1 + 1):length(valid_chan_indices)
            idx2 = valid_chan_indices(i2);
            reg2 = chan_regions{idx2};
            
            if ~strcmp(string(reg1), string(reg2))
                pair_c1 = [pair_c1; idx1]; %#ok<AGROW>
                pair_c2 = [pair_c2; idx2]; %#ok<AGROW>
            end
        end
    end
    
    npairs = length(pair_c1);
    if npairs == 0
        warning('ephys_coherence:NoPairs', 'No inter-regional pairs found for subject %s', sub_name);
        continue;
    end
    
    %% --- 2G. Setup Trial Sorting Conditions ---
    sort_cond_var = op.sort_cond;
    if ~ismember(sort_cond_var, trials_tbl.Properties.VariableNames)
        error('ephys_coherence:CondNotFound', 'Condition variable ''%s'' not in trials table.', sort_cond_var);
    end
    
    if ~isempty(op.sort_cond_vals)
        cond_values = op.sort_cond_vals;
    else
        raw_conds = trials_tbl.(sort_cond_var);
        if iscell(raw_conds) || isstring(raw_conds)
            cond_values = unique(raw_conds(~ismissing(string(raw_conds))), 'stable');
            if isstring(cond_values)
                cond_values = cellstr(cond_values);
            end
        else
            cond_values = unique(raw_conds(~isnan(raw_conds)), 'stable');
            cond_values = num2cell(cond_values);
        end
    end
    
    nconds = length(cond_values);
    
    %% --- 2H. Setup Synchronization Events ---
    sync_tbl = op.sync_events;
    if iscell(sync_tbl.event) || isstring(sync_tbl.event)
        event_names = cellstr(sync_tbl.event);
    else
        event_names = cellstr(sync_tbl.Properties.RowNames);
    end
    nevents = length(event_names);
    
    %% --- 2I. Prepare Output Metadata Columns for Subject Pairs ---
    sub_col = cell(npairs, 1);
    [sub_col{:}] = deal(sub_name);
    
    % Expose chan and region directly as Npairs x 2 string matrices so MATLAB displays values in top-level table view
    chan_col = strings(npairs, 2);
    region_col = strings(npairs, 2);
    
    % n_good_trials stored as Npairs x nconds matrix (trials per condition usable jointly by both electrodes)
    n_good_trials_col = zeros(npairs, nconds);
    
    % Initialize copied resp variables (1 x 2 values per pair)
    resp_copies = struct();
    for ivar = 1:length(op.resp_vars_to_copy)
        vname = op.resp_vars_to_copy{ivar};
        resp_copies.(vname) = zeros(npairs, 2);
    end
    
    for ipair = 1:npairs
        idx1 = pair_c1(ipair);
        idx2 = pair_c2(ipair);
        
        name1 = D.label{idx1};
        name2 = D.label{idx2};
        
        elec1_label = '';
        elec2_label = '';
        if isfield(D, 'electrode_label') && length(D.electrode_label) >= idx1
            elec1_label = D.electrode_label{idx1};
        end
        if isfield(D, 'electrode_label') && length(D.electrode_label) >= idx2
            elec2_label = D.electrode_label{idx2};
        end
        
        chan_col(ipair, :) = [string(name1), string(name2)];
        region_col(ipair, :) = [string(chan_regions{idx1}), string(chan_regions{idx2})];
        
        % Copy channel metrics from resp table for channel 1 and channel 2
        for ivar = 1:length(op.resp_vars_to_copy)
            vname = op.resp_vars_to_copy{ivar};
            val1 = extract_resp_val(resp_tbl, name1, vname, elec1_label);
            val2 = extract_resp_val(resp_tbl, name2, vname, elec2_label);
            resp_copies.(vname)(ipair, :) = [val1, val2];
        end
    end
    
    %% --- 2J. Pairwise Time-Frequency Coherence Computation ---
    n_cycles = 7;
    target_freqs = op.freqs(:)';
    nfreqs = length(target_freqs);
    
    coh_nested_col = cell(npairs, 1);
    
    for ipair = 1:npairs
        idx1 = pair_c1(ipair);
        idx2 = pair_c2(ipair);
        
        % Table storing sync-event level rows
        event_rows = cell(nevents, 1);
        
        % Pre-calculate usable trials per condition across all sync events for top-level n_good_trials
        for icond = 1:nconds
            cond_val = cond_values{icond};
            if iscell(trials_tbl.(sort_cond_var)) || isstring(trials_tbl.(sort_cond_var))
                tr_indices = find(strcmp(string(trials_tbl.(sort_cond_var)), string(cond_val)));
            else
                tr_indices = find(trials_tbl.(sort_cond_var) == unwrap_val(cond_val));
            end
            
            n_good_cond = 0;
            for itr = 1:length(tr_indices)
                tr_row = tr_indices(itr);
                is_tr_usable = true;
                
                for ievent = 1:nevents
                    ev_name = event_names{ievent};
                    if ~ismember(ev_name, trials_tbl.Properties.VariableNames)
                        is_tr_usable = false; break;
                    end
                    t_ev = trials_tbl.(ev_name)(tr_row);
                    if isnan(t_ev), is_tr_usable = false; break; end
                    
                    if ismember('start', sync_tbl.Properties.VariableNames)
                        t_start_rel = sync_tbl.start(ievent);
                        t_end_rel = sync_tbl.end(ievent);
                    else
                        t_start_rel = sync_tbl{ev_name, 'start'};
                        t_end_rel = sync_tbl{ev_name, 'end'};
                    end
                    
                    if ~is_artifact_free(t_ev + t_start_rel, t_ev + t_end_rel, D, idx1, idx2, artifact_master, time_vec)
                        is_tr_usable = false; break;
                    end
                end
                
                if is_tr_usable
                    n_good_cond = n_good_cond + 1;
                end
            end
            n_good_trials_col(ipair, icond) = n_good_cond;
        end
        
        for ievent = 1:nevents
            ev_name = event_names{ievent};
            
            if ismember('start', sync_tbl.Properties.VariableNames)
                t_start_rel = sync_tbl.start(ievent);
                t_end_rel = sync_tbl.end(ievent);
            else
                t_start_rel = sync_tbl{ev_name, 'start'};
                t_end_rel = sync_tbl{ev_name, 'end'};
            end
            
            dt = 1 / fs;
            rel_time_vec = t_start_rel:dt:t_end_rel;
            ntimes = length(rel_time_vec);
            
            cond_rows = cell(nconds, 1);
            
            for icond = 1:nconds
                cond_val = cond_values{icond};
                
                if iscell(trials_tbl.(sort_cond_var)) || isstring(trials_tbl.(sort_cond_var))
                    trial_indices = find(strcmp(string(trials_tbl.(sort_cond_var)), string(cond_val)));
                else
                    trial_indices = find(trials_tbl.(sort_cond_var) == unwrap_val(cond_val));
                end
                
                good_trials = [];
                for itrial = 1:length(trial_indices)
                    tr_row = trial_indices(itrial);
                    if ~ismember(ev_name, trials_tbl.Properties.VariableNames)
                        continue;
                    end
                    
                    t_ev = trials_tbl.(ev_name)(tr_row);
                    if isnan(t_ev)
                        continue;
                    end
                    
                    t_win_start = t_ev + t_start_rel;
                    t_win_end = t_ev + t_end_rel;
                    
                    if is_artifact_free(t_win_start, t_win_end, D, idx1, idx2, artifact_master, time_vec)
                        good_trials = [good_trials; tr_row]; %#ok<AGROW>
                    end
                end
                
                n_good = length(good_trials);
                
                if n_good > 0
                    coh_tf_matrix = compute_wavelet_coherence(...
                        D, idx1, idx2, trials_tbl, ev_name, good_trials, ...
                        t_start_rel, t_end_rel, op.buffer_window_sec, fs, ...
                        target_freqs, n_cycles, op.coh_measure, time_vec);
                else
                    coh_tf_matrix = NaN(nfreqs, ntimes);
                end
                
                % Rownames of freqs table set to exact numerical string representation without f_ and Hz
                freq_rownames = arrayfun(@(f) num2str(f), target_freqs, 'UniformOutput', false);
                freqs_tbl = table(target_freqs', coh_tf_matrix, ...
                    'VariableNames', {'freq', 'coh'}, ...
                    'RowNames', freq_rownames);
                
                cond_val_unwrapped = unwrap_val(cond_val);
                if isnumeric(cond_val_unwrapped) || islogical(cond_val_unwrapped)
                    cond_str = sprintf('cond_%g', cond_val_unwrapped);
                else
                    cond_str = char(cond_val_unwrapped);
                end
                
                cond_rows{icond} = table({cond_val_unwrapped}, n_good, {freqs_tbl}, ...
                    'VariableNames', {'cond', 'n_good_trials', 'freqs'}, ...
                    'RowNames', {cond_str});
            end
            
            trialconds_tbl = vertcat(cond_rows{:});
            event_rows{ievent} = table({ev_name}, {rel_time_vec}, {trialconds_tbl}, ...
                'VariableNames', {'event', 'time', 'trialconds'}, ...
                'RowNames', {ev_name});
        end
        
        coh_nested_col{ipair} = vertcat(event_rows{:});
    end
    
    %% --- 2K. Assemble Master Pair Table for Subject ---
    sub_pair_tbl = table(sub_col, chan_col, region_col, n_good_trials_col, ...
        'VariableNames', {'sub', 'chan', 'region', 'n_good_trials'});
    
    % Append copied response variables
    for ivar = 1:length(op.resp_vars_to_copy)
        vname = op.resp_vars_to_copy{ivar};
        sub_pair_tbl.(vname) = resp_copies.(vname);
    end
    
    % Append nested coherence table column
    sub_pair_tbl.coh = coh_nested_col;
    master_pair_tables{isub} = sub_pair_tbl;
end

%% =======================================================================
%% 3. CONCATENATE RESULTS ACROSS SUBJECTS & SAVE
%% =======================================================================
coh = vertcat(master_pair_tables{:});

% Save results if savefile path is provided
if isfield(op, 'savefile') && ~isempty(op.savefile)
    save_dir = fileparts(op.savefile);
    if ~isempty(save_dir) && ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    save(op.savefile, 'coh', 'subs', 'op', '-v7.3');
end

end


%% =======================================================================
%% HELPER FUNCTIONS
%% =======================================================================

function val = unwrap_val(val)
% Unwraps cell elements if nested
while iscell(val) && numel(val) == 1
    val = val{1};
end
end


function val = extract_resp_val(resp_tbl, chan_name, var_name, elec_label)
% Helper to safely extract variable value for a given channel from resp table
val = NaN;
if isempty(resp_tbl) || ~ismember(var_name, resp_tbl.Properties.VariableNames)
    return;
end

% 1. Match row by RowNames using chan_name or elec_label
if ismember(chan_name, resp_tbl.Properties.RowNames)
    val = resp_tbl{chan_name, var_name};
    val = unwrap_val(val);
    return;
end

if nargin >= 4 && ~isempty(elec_label) && ismember(elec_label, resp_tbl.Properties.RowNames)
    val = resp_tbl{elec_label, var_name};
    val = unwrap_val(val);
    return;
end

% 2. Match against common channel column names
chan_cols = {'chan', 'channel', 'label', 'name', 'electrode', 'electrode_label', 'chan_name', 'channel_name'};
found_col = '';
for ic = 1:length(chan_cols)
    if ismember(chan_cols{ic}, resp_tbl.Properties.VariableNames)
        found_col = chan_cols{ic};
        break;
    end
end

if ~isempty(found_col)
    col_data = string(resp_tbl.(found_col));
    
    idx = find(strcmpi(col_data, string(chan_name)), 1);
    if isempty(idx) && nargin >= 4 && ~isempty(elec_label)
        idx = find(strcmpi(col_data, string(elec_label)), 1);
    end
    
    % Flexible fallback matching (stripping non-digits e.g. '01' vs '1')
    if isempty(idx)
        clean_chan = regexprep(string(chan_name), '^[^\d]*', '');
        if ~isempty(clean_chan)
            idx = find(strcmpi(regexprep(col_data, '^[^\d]*', ''), clean_chan), 1);
        end
    end
    
    if ~isempty(idx)
        val = resp_tbl{idx, var_name};
        val = unwrap_val(val);
        return;
    end
end
end


function is_clean = is_artifact_free(t_start, t_end, D, idx1, idx2, artifact_tbl, time_vec)
% Checks if analysis window is free of NaNs and manual artifacts for both channels
is_clean = false;

% 1. Bounds check relative to continuous data recording
if t_start < time_vec(1) || t_end > time_vec(end)
    return;
end

sample_start = find(time_vec >= t_start, 1, 'first');
sample_end = find(time_vec <= t_end, 1, 'last');

if isempty(sample_start) || isempty(sample_end) || sample_start >= sample_end
    return;
end

% 2. Check for NaNs within trial analysis window
sig1 = D.trial{1}(idx1, sample_start:sample_end);
sig2 = D.trial{1}(idx2, sample_start:sample_end);

if any(isnan(sig1)) || any(isnan(sig2))
    return;
end

% 3. Check for artifact overlap in master artifact table
if ~isempty(artifact_tbl) && ismember('label', artifact_tbl.Properties.VariableNames)
    chan1_name = D.label{idx1};
    chan2_name = D.label{idx2};
    
    art_mask1 = strcmp(string(artifact_tbl.label), string(chan1_name));
    art_mask2 = strcmp(string(artifact_tbl.label), string(chan2_name));
    
    if any(art_mask1)
        starts1 = artifact_tbl.starts(art_mask1);
        ends1 = artifact_tbl.ends(art_mask1);
        if any(starts1 < t_end & ends1 > t_start)
            return;
        end
    end
    
    if any(art_mask2)
        starts2 = artifact_tbl.starts(art_mask2);
        ends2 = artifact_tbl.ends(art_mask2);
        if any(starts2 < t_end & ends2 > t_start)
            return;
        end
    end
end

is_clean = true;
end


function coh_matrix = compute_wavelet_coherence(D, idx1, idx2, trials_tbl, ...
    ev_name, good_trials, t_start_rel, t_end_rel, buffer_sec, fs, freqs, n_cycles, coh_measure, time_vec)
% Computes trial-averaged wavelet cross-spectrum and coherence across frequencies

nfreqs = length(freqs);
dt = 1 / fs;
rel_time_vec = t_start_rel:dt:t_end_rel;
ntimes = length(rel_time_vec);
ngood = length(good_trials);

% Padded segment lengths
t_pad_start = t_start_rel - buffer_sec;
t_pad_end = t_end_rel + buffer_sec;
rel_pad_time = t_pad_start:dt:t_pad_end;
npad_times = length(rel_pad_time);

crop_start = find(rel_pad_time >= t_start_rel, 1, 'first');
crop_indices = crop_start:(crop_start + ntimes - 1);

S_xy = zeros(nfreqs, ntimes);
S_xx = zeros(nfreqs, ntimes);
S_yy = zeros(nfreqs, ntimes);

for itr = 1:ngood
    tr_row = good_trials(itr);
    t_ev = trials_tbl.(ev_name)(tr_row);
    
    t1 = t_ev + t_pad_start;
    t2 = t_ev + t_pad_end;
    
    s1 = find(time_vec >= t1, 1, 'first');
    s2 = s1 + npad_times - 1;
    
    if s2 > size(D.trial{1}, 2)
        continue;
    end
    
    raw1 = D.trial{1}(idx1, s1:s2);
    raw2 = D.trial{1}(idx2, s1:s2);
    
    W1 = zeros(nfreqs, ntimes);
    W2 = zeros(nfreqs, ntimes);
    
    for ifreq = 1:nfreqs
        f0 = freqs(ifreq);
        
        st = n_cycles / (2 * pi * f0);
        t_wavelet = -3*st : dt : 3*st;
        wavelet = pi^(-0.25) .* exp(1i * 2 * pi * f0 .* t_wavelet) .* exp(-t_wavelet.^2 / (2 * st^2));
        
        conv1 = conv(raw1, wavelet, 'same');
        conv2 = conv(raw2, wavelet, 'same');
        
        W1(ifreq, :) = conv1(crop_indices);
        W2(ifreq, :) = conv2(crop_indices);
    end
    
    S_xy = S_xy + (W1 .* conj(W2));
    S_xx = S_xx + (abs(W1).^2);
    S_yy = S_yy + (abs(W2).^2);
end

S_xy = S_xy ./ ngood;
S_xx = S_xx ./ ngood;
S_yy = S_yy ./ ngood;

if strcmp(coh_measure, 'imag')
    coh_matrix = abs(imag(S_xy)) ./ sqrt(S_xx .* S_yy);
else
    coh_matrix = (abs(S_xy).^2) ./ (S_xx .* S_yy);
end

end



%{
===========================================================================
PROMPT REPRODUCTION
===========================================================================
-write a matlab script ephys_coherence.m for for computing coherence between particular pairs of electrodes, structured as: 
   coh = ephys_coherence(subs, op) 
-this function is intended to be run on a series of subjects described in the subs table
-subs table will have a variable ‘paths’, which is a nsubjects x 1 cell array containing, for each subject [row] a subject-specific struct with fields:
---fieldtrip
---trials
---artifact
---electrodes
---resp
…
-the goal of this function is to do analysis inspired by kinygon et al 2015, in particular fig 1e….. https://doi.org/10.1016/j.neuroscience.2015.07.069
-use descriptive variable names, informative section descriptions, and thorough detailed commenting; reproduce this prompt as comments at the bottom of the function 
…
there are 5 key files/variables used in this experiment for each subject; these are specified by the struct subs.paths{isubject}... (isubject) is the index/row within subs table
-1. fieldtrip - .mat file containing struct D with standard fieldtrip formatting, using a single ‘trial’ containing all timepoints for the experiment
….. time field is a 1-D vector with times for each of the samples in .trial{t}, which has rows for channels, columns for timepoints
…. in this prompt I call the fieldtrip struct D, but it may have a different name; this name is specified op.fieldtrip_var_name; default to ‘D_ref’
……. additionally, this struct has a nonstandard (for fieldtrip) ‘electrode_label’ field; this is the same size and type as .label, and indicates the row of the electrodes table that this channel is linked to 
-2. trials - .tsv table contains key timepoints for various events in the trial, like speech onset and offset, stim onset and offset; also contains trial conditions that may be used for grouping trials
…. also, this argument should be accepted as a table instead of a filename; if provided as a table, it will be structured the same way as the table in the .tsv would be expected to be structured
-3. artifact - filename may be multiple files; each contains timepoints for each channel specifying where we had manually marked artifactual timepoints
… each row of this table specifies a single artifactual window for one channel
…chan names are the same as in fieldtrip ‘labels’ field, times correspond to fieldtrip ‘time’ field
-4. electrodes - .tsv table listing info about electrodes; each channel in the fieldtrip file is linked to row in this table via fieldtrip stuct field ‘electrode_label’
…… we will use the electrode info in this table to find the area and region of each channel/electrode, which will be used to decide whether to analyze the channel and which other channels to pair it with for coherence analysis
-5. resp - .mat file with table varable ‘resp’; resp contains a row for each channel (same names as D; we take tuning data for each channel from here specified by op.resp_vars)
….
the function you create will be called by a script with a name like ephys_coherence_seq.m (the suffix ‘seq’ is specific to the project we are calling this func for)....
…… which will set the following options as fields in an ‘op’ structure provided to ephys_coherence:
-sync_events - this is a table with rows ‘event’,’start’,’end’; ‘event’ will match the name of an event in the trial table which analysis windows will be synced to; ‘start’ and ‘end’ specify the amount of time around this event to cut the analysis window in each trial - e.g. start = -2 and end = 3 indicates to cut 2sec before that event and 3sec after that event
-freqs; a numerical vector specifying the frequencies at which to compute coherence
-regions_to_analyze: cell array of strings e.g. {‘IFG’,’STN’,’Thal’}
…. the intention is that we are computing coherence between regions…..
….. so in this case, for every ifg/stn or ifg/thal or stn/thal pair, we compute coherence (but not e.g. across ifg/ifg pairs) 
….. so we should compute ever inter-regional pair’s coherence, within each of the windows specified by sync_events table
… these ‘regions’ will be variables in the processed electrodes table
-sort_cond: before analyze pairwise coherence, data should be split up by trials based on condition (eg ‘novel’,’learned’,’native’)...
…. condition will be the name of a variable in the trials table, matching this option op.sort_cond
….. so we should end up with timecourses of coherence for each electrode pair for each freq band for each unique value of this trial sorting condition
-resp_vars_to_copy …. this will be cell array of strings (e.g. {‘rspv’,’p_prod’,’bad_elc’}
…. the idea is that for each channel, you will find its row in the resp table…
….  then for each of these vars i, copy the value of var specified for it in resp{CHANNAME,op.resp_vars.resp_vars{i}} to that chan’s entry in coh
…. but in coh, there will be 2 cols in the corresponding variable (one for each channel of the pair)
-savefile: filename/path for saving results
….
-when starting to analyze a subject, first try to load the variable with name in op.fieldtrip_var_name from the fieldtrip file…
….. if that variable isn’t in the file, throw error, and specify what was expected and what variables were found instead
-after loading the fieldtrip file for a subject, load the electrodes .tsv as a table ‘electrodes’. 
-the fieldtrip structure D should have field electrode_label, which maps the channels of the fieldtrip struct to the names in electrodes.name. 
---then use
    [electrodes, regiondef] = define_brain_regions(electrodes) % electrodes is the variable, not the filename
to create electrodes.region. then you should be able to figure out the region of each channel in the fieldtrip struct by mapping D.label to electrodes row, and finding that row of electrodes.region
-load all of the artifact tables and stack them together as the master artifact table for this subject
…if the artifact file is a string rather than a cell array of strings, then it’s just one file to load rather than multiple files to load and combine
….
-when analyzing coherence for an electrode pair, 
…exclude from analysis all trials where any timepoints within the analysis window for either electrode is marked as artifactual in the artifact table 
…(this might vary by event window, e.g. the stim window might be artifactual while the speech window is not)
…also exclude trials where there are any nans for either channel within the analysis window
….
-provide coh as function output
-save into the file op.savefile these variables: coh, subs, op
-coh is a table with 1 row per channel pair, with variables:
----sub - 1 x npairs cell - subject name that this pair belongs to
----chan - 2 x npairs cell array of strings - the names of the 2 channels in the pair, taken from D.label
----region - 2 x npairs cell
----n_good_trials - 1 x npairs int - the number of trials which are usable jointly from both electrodes in the pair (trials must be usable for both) 
----[the other variables copied from resp table, specified by op.resp_vars_to_copy]

CLARIFICATIONS ADDED:
-use morlet wavelet, with 7 wavelets; comment that kingyon 2015 used short time FFT
-include a param op.coh_measure…. user can choose between ‘mag_sq’ for magnitude squared (default) and ‘imag’ (imaginary)... comment that kingyon 2015 use mag_sq
-include a param op.buffer_window_sec…. default to 1sec
-includ a param op.sort_cond_vals: if provided, this is the list of values of sort_cond to use (ignore others); always store them in this specified order; if not provided, use all vals of sort_cond

how to store coherence analysis in the table variable coh.coh{chan_pair_index}
-coh.coh contains one cell per chan pair (row of coh)
-coh.coh cell contains a table with one row per sync event; thus each row contains the analysis for this pair for e.g. speech-onset-aligned coherence
….. has column ‘time’ which is a cell containing the times relative to sync-time onset (which is x=0) of the datapoints stored deeper; this will be used for plotting as x values while coherence is y or z values
….. has column event with the event label (also use these are rownames); also has column ‘trialconds’ which is cells - going a level deeper to contain trial conditions
-trialconds is a table with row for each value of op.sort_cond (make cond both a table column and the row names)
…. trialconds has a column ‘freqs’ which is a cell containing another level of nested table
-freqs table has a row for each freq listed in op.freqs - make these a table variable and the table rownames
…. freqs also has a variable ‘coh’; coh is a n_freqs x ntimes double, containing the actual coherence values for each timepoint for each freq

-line 45 of ephys_coherence_seq.m is updated to :  coh = ephys_coherence(subs, op)
-”ephys_coherence_seq.m, op.freqs = [20, 100]” …. specifies two discrete frequencies (20hz and 100hz)
===========================================================================
%}