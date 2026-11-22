function [sync_out, op_out] = sync_audio_and_annot(sync_in, op)
% SYNC_AUDIO_AND_ANNOT Align moving audio/video recordings to a master audio
% clock, update annotation tables with master-synced timestamps, and export
% buffered event-wise audio clips.
%
% USAGE:
%   [sync_out, op_out] = sync_audio_and_annot(sync_in, op)
%
% INPUTS:
%   sync_in : Path to a sync-landmarks table (.tsv or .xlsx) OR a MATLAB
%             table variable containing master and moving recording rows.
%   op      : Struct of configuration options (all fields optional; defaults
%             applied automatically and returned in op_out).
%
% =========================================================================
% OPTIONS STRUCT (op) PARAMETERS AND DEFAULTS
% =========================================================================
%   --- Grouping & File Filtering ---
%   op.grouping_columns         : (default: 'step') Column name (char/string)
%                                 or cell array of column names (e.g.,
%                                 {'step','task'} or {'task','run'}) defining
%                                 recording time blocks. Each group must have
%                                 exactly 1 row with master == 1.
%   op.include_criteria         : (default: {}) Cell array of even length
%                                 {'col1','val1','col2','val2',...} specifying
%                                 which moving rows to include. Empty = all.
%   op.exclude_criteria         : (default: {}) Cell array of even length
%                                 {'col1','val1',...} specifying moving rows
%                                 to exclude (takes precedence over include).
%
%   --- Paths & Output Naming ---
%   op.sync_table_savename      : (default: '') Output path (.tsv or .xlsx)
%                                 for synced table. If empty, appends '_sync'
%                                 to sync_in filename, or saves to
%                                 ./sync_<timestamp>.xlsx if sync_in is a table.
%   op.annot_dir                : (default: {}) Char/string path or cell array
%                                 of paths where annotation files live. Also
%                                 searches one level down in 'unsynced/'. If
%                                 empty, uses moving audio directories.
%   op.annotation_file_suffixes : (default: {''}) Cell array of suffixes
%                                 (e.g., {'scored','annotated',''}) matched on
%                                 annotation filenames (with optional '_' or
%                                 '-' prefix before the suffix).
%   op.eventwise_audvid_dir     : (default: '') Top-level directory for chopped
%                                 event-wise audio clips (also accepts alias
%                                 op.event_wise_audvid_dir). If empty, audio
%                                 chopping is skipped.
%   op.event_audvid_buffer_sec  : (default: [1 1]) 2-element vector [pre post]
%                                 in seconds specifying buffer around each
%                                 annotated event when chopping audio.
%
%   --- Signal Conditioning & Two-Stage Alignment ---
%   op.target_fs                : (default: 10000) Common sample rate (Hz)
%                                 used for cross-correlation.
%   op.bandpass_hz              : (default: [100 4500]) Zero-phase bandpass
%                                 filter cutoffs [low high] in Hz.
%   op.env_win_ms               : (default: 25) RMS envelope window (ms) for
%                                 Stage 1 coarse macro-alignment.
%   op.env_fs                   : (default: 200) Envelope sample rate (Hz) for
%                                 Stage 1 macro-alignment.
%   op.local_z_sec              : (default: 10) Sliding z-score window (s) to
%                                 equalize mic level differences in Stage 1.
%   op.fine_win_sec             : (default: 20) Window length (s) for Stage 2
%                                 piecewise local cross-correlation.
%   op.fine_step_sec            : (default: 30) Step interval (s) between
%                                 Stage 2 local correlation windows.
%   op.fine_env_smooth_ms       : (default: 2) Rectified signal smoothing (ms)
%                                 for Stage 2 fine correlation (handles room
%                                 mic acoustic phase differences).
%   op.fine_search_rad_seeded_sec   : (default: 0.5) Local lag search radius
%                                 (s) around manual t1/t2 seed trajectory.
%   op.fine_search_rad_unseeded_sec : (default: 1.5) Local lag search radius
%                                 (s) around Stage 1 global lag when unseeded.
%   op.inlier_tol_ms            : (default: 20) Residual threshold (ms) for
%                                 counting a Stage 2 window as an inlier.
%
%   --- Quality Control (QC) Cutoffs & Plotting ---
%   op.qc_min_corr              : (default: 0.25) Minimum median active-window
%                                 Pearson correlation (r_median).
%   op.qc_min_peak_z            : (default: 6.0) Minimum peak-to-sidelobe
%                                 Z-score (Z_peak).
%   op.qc_min_inlier_frac       : (default: 0.50) Minimum fraction of active
%                                 windows agreeing with linear drift fit.
%   op.qc_min_inlier_win        : (default: 4) Minimum count of inlier windows.
%   op.qc_max_resid_ms          : (default: 15.0) Maximum inlier residual
%                                 RMSE (ms) from linear clock model.
%   op.qc_max_warp_ppm          : (default: 2000) Maximum clock drift magnitude
%                                 |warpfactor - 1| * 1e6 (ppm).
%   op.qc_max_landmark_err_ms   : (default: 250) Maximum discrepancy (ms)
%                                 between auto-alignment and manual t1/t2
%                                 landmarks before prompting the user.
%   op.qc_save_figs             : (default: true) Save 3-panel diagnostic QC
%                                 PNG figure for each aligned pair.
%   op.show_qc_figs             : (default: false) Display QC figures on screen
%                                 (true/'on' or false/'off').
%   op.qc_figs_savedir          : (default: '') Directory to save QC figures.
%                                 If empty, creates 'sync_qc_figs' inside the
%                                 input sync_in file's folder (or pwd if
%                                 sync_in was passed as a table variable).
%
% OUTPUTS:
%   sync_out : Updated sync table containing 'warpfactor', 'offset', and QC
%              metric columns ('sync_status', 'r_median', 'clock_drift_ppm',
%              'landmark_err_t1_ms', 'landmark_err_t2_ms').
%   op_out   : Options struct populated with all defaults used.

%% 1. Parse Options & Load Sync Table
if nargin < 2 || isempty(op)
    op = struct();
end
[sync_tbl, sync_source_dir, sync_source_file] = load_sync_table(sync_in);
op_out = parse_and_validate_options(op, sync_tbl, sync_source_dir, sync_source_file);

n_rows = height(sync_tbl);

% Initialize output columns on sync_out
sync_out = sync_tbl;
sync_out.warpfactor         = nan(n_rows, 1);
sync_out.offset             = nan(n_rows, 1);
sync_out.sync_status        = repmat({'unprocessed'}, n_rows, 1);
sync_out.r_median           = nan(n_rows, 1);
sync_out.clock_drift_ppm    = nan(n_rows, 1);
sync_out.landmark_err_t1_ms = nan(n_rows, 1);
sync_out.landmark_err_t2_ms = nan(n_rows, 1);

% Extract master flag and landmark columns safely
master_col = get_numeric_col(sync_tbl, 'master');
has_t_cols = all(ismember({'t1', 't2'}, sync_tbl.Properties.VariableNames));
if has_t_cols
    t1_col = get_numeric_col(sync_tbl, 't1');
    t2_col = get_numeric_col(sync_tbl, 't2');
else
    t1_col = nan(n_rows, 1);
    t2_col = nan(n_rows, 1);
end

% Tracking structs for end-of-run summary
aligned_log = struct('group', {}, 'moving_file', {}, 'master_file', {}, ...
    'method', {}, 'warpfactor', {}, 'offset', {}, 'r_median', {}, ...
    'clock_drift_ppm', {}, 'err_t1_ms', {}, 'err_t2_ms', {}, ...
    'qc_passed', {}, 'annot_files', {}, 'n_clips', {});
skipped_log = struct('group', {}, 'file', {}, 'reason', {});

%% 2. Group Rows & Process Each Recording Block
[group_ids, group_labels] = resolve_group_keys(sync_tbl, op_out.grouping_columns);
n_groups = numel(group_labels);

fprintf('\n====================================================================\n');
fprintf('Starting audio & annotation sync across %d group(s)...\n', n_groups);
fprintf('====================================================================\n');

for g = 1:n_groups
    g_rows = find(group_ids == g);
    g_label = group_labels{g};
    is_master = (master_col(g_rows) == 1);
    n_masters = sum(is_master);

    % Validate master count per group
    if n_masters == 0
        msg = sprintf('Group [%s] has 0 master rows (master==1); skipping group.', g_label);
        warning('sync_audio_and_annot:NoMaster', '%s', msg);
        for idx = g_rows(:)'
            fname = get_display_filename(sync_tbl, idx);
            sync_out.sync_status{idx} = 'skipped_no_master';
            skipped_log(end+1) = struct('group', g_label, 'file', fname, ...
                'reason', 'No master file (master==1) in group'); %#ok<AGROW>
        end
        continue;
    elseif n_masters > 1
        error('sync_audio_and_annot:MultipleMasters', ...
            'Group [%s] has %d master rows (expected exactly 1).', g_label, n_masters);
    end

    master_idx = g_rows(is_master);
    moving_indices = g_rows(~is_master);

    % Set master clock parameters
    sync_out.warpfactor(master_idx)  = 1;
    sync_out.offset(master_idx)      = 0;
    sync_out.sync_status{master_idx} = 'master';

    if isempty(moving_indices)
        continue;
    end

    % Resolve master audio file
    master_path = resolve_media_path(sync_tbl, master_idx, sync_source_dir);
    master_disp = get_display_filename(sync_tbl, master_idx);
    if isempty(master_path) || ~isfile(master_path)
        reason = sprintf('Master audio file not found: %s', master_path);
        warning('sync_audio_and_annot:MissingMasterFile', '%s', reason);
        sync_out.sync_status{master_idx} = 'master_file_missing';
        for idx = moving_indices(:)'
            sync_out.sync_status{idx} = 'skipped_missing_master_file';
            skipped_log(end+1) = struct('group', g_label, ...
                'file', get_display_filename(sync_tbl, idx), 'reason', reason); %#ok<AGROW>
        end
        continue;
    end

    % Lazy-load conditioned master audio only if at least one moving row passes filters
    master_audio_loaded = false;
    y_mas = [];
    mas_dur = NaN;

    for m_idx = moving_indices(:)'
        mov_disp = get_display_filename(sync_tbl, m_idx);

        % Check include / exclude criteria
        [keep_row, filter_reason] = check_criteria_match(sync_tbl, m_idx, ...
            op_out.include_criteria, op_out.exclude_criteria);
        if ~keep_row
            sync_out.sync_status{m_idx} = 'skipped_filtered';
            skipped_log(end+1) = struct('group', g_label, 'file', mov_disp, ...
                'reason', filter_reason); %#ok<AGROW>
            continue;
        end

        % Resolve moving audio file (including video -> .wav lookup + duration check)
        [mov_audio_path, orig_media_path, resolve_ok, resolve_reason] = ...
            resolve_moving_audio_path(sync_tbl, m_idx, sync_source_dir);
        if ~resolve_ok
            sync_out.sync_status{m_idx} = 'skipped_media_error';
            skipped_log(end+1) = struct('group', g_label, 'file', mov_disp, ...
                'reason', resolve_reason); %#ok<AGROW>
            continue;
        end

        % Load master audio on first valid moving pair in this group
        if ~master_audio_loaded
            fprintf('\n[%s] Loading master audio: %s\n', g_label, master_disp);
            [y_mas, mas_dur] = load_and_condition_audio(master_path, op_out);
            master_audio_loaded = true;
        end

        fprintf('[%s] Aligning moving audio: %s -> %s\n', g_label, mov_disp, master_disp);
        [y_mov, mov_dur] = load_and_condition_audio(mov_audio_path, op_out);

        % Check for valid manual landmarks on both master and moving rows
        landmarks = struct('valid', false, 't_mas_1', NaN, 't_mas_2', NaN, ...
            't_mov_1', NaN, 't_mov_2', NaN, 'w_manual', NaN, 'off_manual', NaN);
        if has_t_cols
            tm1 = t1_col(master_idx); tm2 = t2_col(master_idx);
            tv1 = t1_col(m_idx);      tv2 = t2_col(m_idx);
            if all(isfinite([tm1, tm2, tv1, tv2])) && (tm2 > tm1) && (tv2 > tv1)
                if tm1 >= 0 && tm2 <= (mas_dur + 0.5) && tv1 >= 0 && tv2 <= (mov_dur + 0.5)
                    landmarks.valid      = true;
                    landmarks.t_mas_1    = tm1;
                    landmarks.t_mas_2    = tm2;
                    landmarks.t_mov_1    = tv1;
                    landmarks.t_mov_2    = tv2;
                    landmarks.w_manual   = (tm2 - tm1) / (tv2 - tv1);
                    landmarks.off_manual = tm1 - landmarks.w_manual * tv1;
                else
                    warning('sync_audio_and_annot:LandmarkOutOfBounds', ...
                        'Landmarks for %s exceed file durations; running unseeded.', mov_disp);
                end
            end
        end

        % Run two-stage audio alignment
        res = align_audio_pair(y_mas, y_mov, op_out, landmarks);

        % Save/display diagnostic QC figure
        if op_out.qc_save_figs || is_true_flag(op_out.show_qc_figs)
            if ~isfolder(op_out.qc_figs_savedir)
                mkdir(op_out.qc_figs_savedir);
            end
            [~, mov_stem] = fileparts(mov_audio_path);
            fig_path = fullfile(op_out.qc_figs_savedir, [mov_stem, '_syncQC.png']);
            plot_alignment_qc(y_mas, y_mov, res, landmarks, op_out, ...
                g_label, mov_disp, master_disp, fig_path);
        end

        % Evaluate QC decision & handle >250 ms manual landmark discrepancy
        align_success = false;
        align_method  = 'auto';
        final_w       = NaN;
        final_off     = NaN;
        fail_reason   = '';

        if ~res.audio_qc_passed
            fail_reason = sprintf('Failed audio QC (%s)', strjoin(res.qc_fail_reasons, '; '));
            sync_out.sync_status{m_idx} = 'failed_qc';
        elseif landmarks.valid && (res.landmark_max_err_ms > op_out.qc_max_landmark_err_ms)
            % Prompt user when manual landmarks disagree with auto-alignment by >250 ms
            fprintf('\n*** MANUAL LANDMARK DISCREPANCY ALERT [%s] ***\n', g_label);
            fprintf('  Moving file : %s\n', mov_disp);
            fprintf('  Master file : %s\n', master_disp);
            fprintf('  Auto-align  : warpfactor = %.8f (%.1f ppm), offset = %.4f s, r_median = %.3f\n', ...
                res.warpfactor, res.clock_drift_ppm, res.offset, res.r_median);
            fprintf('  Manual t1/t2: warpfactor = %.8f (%.1f ppm), offset = %.4f s\n', ...
                landmarks.w_manual, (landmarks.w_manual - 1)*1e6, landmarks.off_manual);
            fprintf('  Endpoint discrepancy: dt1 = %+.1f ms, dt2 = %+.1f ms (threshold: %.1f ms)\n', ...
                res.landmark_err_t1_ms, res.landmark_err_t2_ms, op_out.qc_max_landmark_err_ms);

            choice = '';
            while ~ismember(choice, {'m', 'f'})
                raw_in = input('  Choose [m] to fall back to manual t1/t2 solution, or [f] to mark sync as FAILED: ', 's');
                choice = lower(strtrim(raw_in));
                if strcmp(choice, 'manual'), choice = 'm'; end
                if ismember(choice, {'fail', 'failed', 'reject'}), choice = 'f'; end
            end

            if strcmp(choice, 'm')
                align_success = true;
                align_method  = 'manual_fallback';
                final_w       = landmarks.w_manual;
                final_off     = landmarks.off_manual;
                sync_out.sync_status{m_idx} = 'aligned_manual_fallback';
            else
                align_success = false;
                fail_reason   = sprintf('User rejected landmark mismatch (dt1=%+.1f ms, dt2=%+.1f ms)', ...
                    res.landmark_err_t1_ms, res.landmark_err_t2_ms);
                sync_out.sync_status{m_idx} = 'failed_landmark_mismatch';
            end
        else
            align_success = true;
            align_method  = 'auto';
            final_w       = res.warpfactor;
            final_off     = res.offset;
            sync_out.sync_status{m_idx} = 'aligned_auto';
        end

        % Record QC metrics in sync_out (save r_median, clock drift, and endpoint errors even if rejected)
        sync_out.r_median(m_idx)           = res.r_median;
        sync_out.landmark_err_t1_ms(m_idx) = res.landmark_err_t1_ms;
        sync_out.landmark_err_t2_ms(m_idx) = res.landmark_err_t2_ms;

        if align_success
            final_drift_ppm = (final_w - 1) * 1e6;
            sync_out.warpfactor(m_idx)      = final_w;
            sync_out.offset(m_idx)          = final_off;
            sync_out.clock_drift_ppm(m_idx) = final_drift_ppm;

            % Display QC metrics in console
            fprintf('  -> SUCCESS (%s): warpfactor=%.8f, offset=%.4f s | r_median=%.3f, drift=%+.1f ppm', ...
                align_method, final_w, final_off, res.r_median, final_drift_ppm);
            if landmarks.valid
                fprintf(', dt1=%+.1f ms, dt2=%+.1f ms\n', res.landmark_err_t1_ms, res.landmark_err_t2_ms);
            else
                fprintf(' (no manual landmarks)\n');
            end

            % Process associated annotation files & chop event-wise audio
            [synced_annots, total_clips] = process_annotation_files( ...
                mov_audio_path, orig_media_path, final_w, final_off, op_out);

            aligned_log(end+1) = struct( ...
                'group', g_label, ...
                'moving_file', mov_disp, ...
                'master_file', master_disp, ...
                'method', align_method, ...
                'warpfactor', final_w, ...
                'offset', final_off, ...
                'r_median', res.r_median, ...
                'clock_drift_ppm', final_drift_ppm, ...
                'err_t1_ms', res.landmark_err_t1_ms, ...
                'err_t2_ms', res.landmark_err_t2_ms, ...
                'qc_passed', true, ...
                'annot_files', {synced_annots}, ...
                'n_clips', total_clips); %#ok<AGROW>
        else
            sync_out.clock_drift_ppm(m_idx) = res.clock_drift_ppm;
            fprintf('  -> FAILED: %s | r_median=%.3f, est_drift=%+.1f ppm', ...
                fail_reason, res.r_median, res.clock_drift_ppm);
            if landmarks.valid
                fprintf(', dt1=%+.1f ms, dt2=%+.1f ms\n', res.landmark_err_t1_ms, res.landmark_err_t2_ms);
            else
                fprintf('\n');
            end
            skipped_log(end+1) = struct('group', g_label, 'file', mov_disp, ...
                'reason', fail_reason); %#ok<AGROW>
        end
    end
end

%% 3. Save Updated Sync Table & Print End-of-Run Summary
save_dir = fileparts(op_out.sync_table_savename);
if ~isempty(save_dir) && ~isfolder(save_dir)
    mkdir(save_dir);
end
[~, ~, save_ext] = fileparts(op_out.sync_table_savename);
if strcmpi(save_ext, '.tsv') || strcmpi(save_ext, '.txt')
    writetable(sync_out, op_out.sync_table_savename, 'FileType', 'text', 'Delimiter', '\t');
else
    writetable(sync_out, op_out.sync_table_savename);
end

print_sync_summary(aligned_log, skipped_log, op_out);

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function op = parse_and_validate_options(op, sync_tbl, sync_source_dir, sync_source_file)
% Set default values for all options and validate filter / path arguments.

% Support alias op.event_wise_audvid_dir -> op.eventwise_audvid_dir
if isfield(op, 'event_wise_audvid_dir') && (~isfield(op, 'eventwise_audvid_dir') || isempty(op.eventwise_audvid_dir))
    op.eventwise_audvid_dir = op.event_wise_audvid_dir;
end

% Default grouping_columns: 'step' if present, else {'task','run'} if present
if ~isfield(op, 'grouping_columns') || isempty(op.grouping_columns)
    if ismember('step', sync_tbl.Properties.VariableNames)
        op.grouping_columns = 'step';
    elseif all(ismember({'task', 'run'}, sync_tbl.Properties.VariableNames))
        op.grouping_columns = {'task', 'run'};
    else
        op.grouping_columns = 'step';
    end
end
if ischar(op.grouping_columns) || isstring(op.grouping_columns)
    op.grouping_columns = cellstr(op.grouping_columns);
end

% Filtering criteria
if ~isfield(op, 'include_criteria'), op.include_criteria = {}; end
if ~isfield(op, 'exclude_criteria'), op.exclude_criteria = {}; end

if isstring(op.include_criteria), op.include_criteria = cellstr(op.include_criteria); end
if isstring(op.exclude_criteria), op.exclude_criteria = cellstr(op.exclude_criteria); end

if mod(numel(op.include_criteria), 2) ~= 0
    error('sync_audio_and_annot:InvalidCriteria', 'op.include_criteria must have an even number of elements (col, val pairs).');
end
if mod(numel(op.exclude_criteria), 2) ~= 0
    error('sync_audio_and_annot:InvalidCriteria', 'op.exclude_criteria must have an even number of elements (col, val pairs).');
end
if ~isempty(op.include_criteria) && ~isempty(op.exclude_criteria) && isequal(op.include_criteria, op.exclude_criteria)
    error('sync_audio_and_annot:IdenticalCriteria', ...
        'op.include_criteria and op.exclude_criteria are non-empty and identical.');
end

% Paths & annotation options
if ~isfield(op, 'annotation_file_suffixes') || isempty(op.annotation_file_suffixes)
    op.annotation_file_suffixes = {''};
elseif ischar(op.annotation_file_suffixes) || isstring(op.annotation_file_suffixes)
    op.annotation_file_suffixes = cellstr(op.annotation_file_suffixes);
end

if ~isfield(op, 'annot_dir') || isempty(op.annot_dir)
    op.annot_dir = {};
elseif ischar(op.annot_dir) || isstring(op.annot_dir)
    op.annot_dir = cellstr(op.annot_dir);
end

if ~isfield(op, 'sync_table_savename') || isempty(op.sync_table_savename)
    if ~isempty(sync_source_file)
        [p, n, ext] = fileparts(sync_source_file);
        op.sync_table_savename = fullfile(p, [n, '_sync', ext]);
    else
        stamp = datestr(now, 'yyyymmdd_HHMMSS');
        op.sync_table_savename = fullfile(pwd, sprintf('sync_%s.xlsx', stamp));
    end
end

if ~isfield(op, 'eventwise_audvid_dir'), op.eventwise_audvid_dir = ''; end
if ~isfield(op, 'event_audvid_buffer_sec') || isempty(op.event_audvid_buffer_sec)
    op.event_audvid_buffer_sec = [1 1];
end

% Alignment algorithm parameters
defaults = struct( ...
    'target_fs', 10000, ...
    'bandpass_hz', [100 4500], ...
    'env_win_ms', 25, ...
    'env_fs', 200, ...
    'local_z_sec', 10, ...
    'fine_win_sec', 20, ...
    'fine_step_sec', 30, ...
    'fine_env_smooth_ms', 2, ...
    'fine_search_rad_seeded_sec', 0.5, ...
    'fine_search_rad_unseeded_sec', 1.5, ...
    'inlier_tol_ms', 20, ...
    'qc_min_corr', 0.25, ...
    'qc_min_peak_z', 6.0, ...
    'qc_min_inlier_frac', 0.50, ...
    'qc_min_inlier_win', 4, ...
    'qc_max_resid_ms', 15.0, ...
    'qc_max_warp_ppm', 2000, ...
    'qc_max_landmark_err_ms', 250, ...
    'qc_save_figs', true, ...
    'show_qc_figs', false, ...
    'qc_figs_savedir', '');

fnames = fieldnames(defaults);
for k = 1:numel(fnames)
    fn = fnames{k};
    if ~isfield(op, fn) || isempty(op.(fn))
        op.(fn) = defaults.(fn);
    end
end

% Default qc_figs_savedir: 'sync_qc_figs' next to sync_in file (or in pwd if table input)
if isempty(op.qc_figs_savedir)
    if ~isempty(sync_source_dir)
        op.qc_figs_savedir = fullfile(sync_source_dir, 'sync_qc_figs');
    else
        op.qc_figs_savedir = fullfile(pwd, 'sync_qc_figs');
    end
end
end

% -------------------------------------------------------------------------
function [tbl, source_dir, source_file] = load_sync_table(sync_in)
if istable(sync_in)
    tbl = sync_in;
    source_dir  = '';
    source_file = '';
elseif ischar(sync_in) || isstring(sync_in)
    source_file = char(sync_in);
    if ~isfile(source_file)
        error('sync_audio_and_annot:FileNotFound', 'Sync table file not found: %s', source_file);
    end
    source_dir = fileparts(source_file);
    [~, ~, ext] = fileparts(source_file);
    if strcmpi(ext, '.tsv') || strcmpi(ext, '.txt')
        tbl = readtable(source_file, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
    elseif strcmpi(ext, '.xlsx') || strcmpi(ext, '.xls')
        tbl = readtable(source_file, 'TextType', 'string');
    else
        error('sync_audio_and_annot:UnsupportedFormat', ...
            'Unsupported sync table extension "%s" (expected .tsv or .xlsx).', ext);
    end
else
    error('sync_audio_and_annot:InvalidInput', 'sync_in must be a file path (.tsv/.xlsx) or a MATLAB table.');
end

if ~ismember('master', tbl.Properties.VariableNames)
    error('sync_audio_and_annot:MissingColumn', 'Sync table must contain a "master" column.');
end
end

% -------------------------------------------------------------------------
function [group_ids, group_labels] = resolve_group_keys(tbl, grouping_cols)
missing = setdiff(grouping_cols, tbl.Properties.VariableNames);
if ~isempty(missing)
    error('sync_audio_and_annot:MissingGroupCol', ...
        'Grouping column(s) not found in sync table: %s', strjoin(missing, ', '));
end

n = height(tbl);
key_strs = strings(n, 1);
for c = 1:numel(grouping_cols)
    col_name = grouping_cols{c};
    vals = string(get_str_col(tbl, col_name));
    part = string(col_name) + "=" + vals;
    if c == 1
        key_strs = part;
    else
        key_strs = key_strs + " | " + part;
    end
end

[unique_keys, ~, group_ids] = unique(key_strs, 'stable');
group_labels = cellstr(unique_keys);
end

% -------------------------------------------------------------------------
function [keep, reason] = check_criteria_match(tbl, row_idx, inc_crit, exc_crit)
keep = true;
reason = '';

% 1. Check include_criteria (AND across distinct columns, OR within repeated columns)
if ~isempty(inc_crit)
    inc_cols = inc_crit(1:2:end);
    inc_vals = inc_crit(2:2:end);
    u_cols = unique(inc_cols, 'stable');
    for c = 1:numel(u_cols)
        cname = u_cols{c};
        if ~ismember(cname, tbl.Properties.VariableNames)
            keep = false;
            reason = sprintf('Include column "%s" missing from table', cname);
            return;
        end
        row_val = strtrim(get_cell_str_at(tbl, cname, row_idx));
        allowed = inc_vals(strcmp(inc_cols, cname));
        matched = any(cellfun(@(v) strcmpi(row_val, strtrim(char(string(v)))), allowed));
        if ~matched
            keep = false;
            reason = sprintf('Did not match include_criteria (%s = %s)', cname, row_val);
            return;
        end
    end
end

% 2. Check exclude_criteria (takes precedence; any match excludes row)
if ~isempty(exc_crit)
    exc_cols = exc_crit(1:2:end);
    exc_vals = exc_crit(2:2:end);
    for k = 1:numel(exc_cols)
        cname = exc_cols{k};
        if ismember(cname, tbl.Properties.VariableNames)
            row_val = strtrim(get_cell_str_at(tbl, cname, row_idx));
            target_val = strtrim(char(string(exc_vals{k})));
            if strcmpi(row_val, target_val)
                keep = false;
                reason = sprintf('Matched exclude_criteria (%s = %s)', cname, row_val);
                return;
            end
        end
    end
end
end

% -------------------------------------------------------------------------
function full_p = resolve_media_path(tbl, row_idx, base_dir)
vnames = tbl.Properties.VariableNames;
fp = '';
fn = '';

dir_cols = {'dir', 'filepath', 'fullpath', 'path', 'folder', 'directory'};
for c = 1:numel(dir_cols)
    if ismember(dir_cols{c}, vnames)
        val = strtrim(get_cell_str_at(tbl, dir_cols{c}, row_idx));
        if ~isempty(val) && ~strcmpi(val, 'NaN')
            fp = val;
            break;
        end
    end
end

file_cols = {'filename', 'file', 'name'};
for c = 1:numel(file_cols)
    if ismember(file_cols{c}, vnames)
        val = strtrim(get_cell_str_at(tbl, file_cols{c}, row_idx));
        if ~isempty(val) && ~strcmpi(val, 'NaN')
            fn = val;
            break;
        end
    end
end

if ~isempty(fp) && ~isempty(fn)
    [~, ~, fp_ext] = fileparts(fp);
    if isempty(fp_ext)
        full_p = fullfile(fp, fn);
    else
        full_p = fp;
    end
elseif ~isempty(fp)
    full_p = fp;
else
    full_p = fn;
end

if isempty(full_p) || isfile(full_p)
    return;
end

% Normalize path separators
full_p_norm = strrep(strrep(full_p, '\', filesep), '/', filesep);
if isfile(full_p_norm)
    full_p = full_p_norm;
    return;
end

if ~isempty(base_dir)
    % 1. Check relative to sync_in directory and up to 4 ancestor levels (e.g., R01-SML_data_shared)
    anc = base_dir;
    for up = 0:4
        if up > 0
            parent_anc = fileparts(anc);
            if isempty(parent_anc) || strcmp(parent_anc, anc), break; end
            anc = parent_anc;
        end
        cand = fullfile(anc, full_p_norm);
        if isfile(cand)
            full_p = cand;
            return;
        end
    end

    % 2. Remap shared dataset root across machines (e.g. C:\Dropbox\R01-SML_data_shared -> local Dropbox path)
    markers = {'R01-SML_data_shared', 'sourcedata', 'derivatives'};
    for m = 1:numel(markers)
        mk = markers{m};
        idx_p = strfind(lower(full_p_norm), lower(mk));
        idx_b = strfind(lower(base_dir), lower(mk));
        if ~isempty(idx_p) && ~isempty(idx_b)
            remapped = [base_dir(1:idx_b(1)-1), full_p_norm(idx_p(1):end)];
            if isfile(remapped)
                full_p = remapped;
                return;
            end
        end
    end
end
end









% -------------------------------------------------------------------------
function [audio_path, orig_path, ok, reason] = resolve_moving_audio_path(tbl, row_idx, base_dir)
orig_path = resolve_media_path(tbl, row_idx, base_dir);
audio_path = orig_path;
ok = false;
reason = '';

if isempty(orig_path)
    reason = 'Empty filepath/filename in sync table';
    return;
end

[~, ~, ext] = fileparts(orig_path);
video_exts = {'.mp4', '.mov', '.avi', '.mkv', '.m4v', '.mts', '.mpg', '.mpeg', '.wmv'};
is_video = any(strcmpi(ext, video_exts));
if ~is_video && ismember('filetype', tbl.Properties.VariableNames)
    ft = lower(strtrim(get_cell_str_at(tbl, 'filetype', row_idx)));
    if contains(ft, 'video') || any(contains(ft, {'mp4', 'mov', 'avi', 'mkv'}))
        is_video = true;
    end
end

if is_video
    [p, n, ~] = fileparts(orig_path);
    audio_path = fullfile(p, [n, '.wav']);
    if ~isfile(audio_path)
        reason = sprintf('Extracted .wav audio not found for video: %s', audio_path);
        return;
    end
    % Basic duration consistency check between video and extracted .wav
    try
        ainfo = audioinfo(audio_path);
        wav_dur = ainfo.Duration;
        if isfile(orig_path)
            vobj = VideoReader(orig_path);
            vid_dur = vobj.Duration;
            if abs(vid_dur - wav_dur) > 1.0
                reason = sprintf('Video (%.2fs) and extracted .wav (%.2fs) lengths differ by %.2fs', ...
                    vid_dur, wav_dur, abs(vid_dur - wav_dur));
                return;
            end
        end
    catch ME
        warning('sync_audio_and_annot:VideoDurationCheck', ...
            'Could not verify video duration for %s (%s); proceeding with .wav.', orig_path, ME.message);
    end
else
    if ~isfile(audio_path)
        reason = sprintf('Moving audio file not found: %s', audio_path);
        return;
    end
end

ok = true;
end

% -------------------------------------------------------------------------
function [y_cond, dur_sec] = load_and_condition_audio(filepath, op)
[y_raw, fs_raw] = audioread(filepath);
dur_sec = size(y_raw, 1) / fs_raw;

% Convert stereo/multi-channel to mono (pick cleaner channel if one is dead/clipped)
if size(y_raw, 2) > 1
    rms_ch = rms(y_raw, 1);
    if max(rms_ch) > 10 * min(rms_ch) && min(rms_ch) < 1e-4
        [~, best_ch] = max(rms_ch);
        y_mono = y_raw(:, best_ch);
    else
        y_mono = mean(y_raw, 2);
    end
else
    y_mono = y_raw(:);
end
clear y_raw;

% Resample to common target_fs
if fs_raw ~= op.target_fs
    [p, q] = rat(op.target_fs / fs_raw, 1e-6);
    if exist('resample', 'file')
        y_rs = resample(double(y_mono), p, q);
    else
        t_old = (0:numel(y_mono)-1)' / fs_raw;
        t_new = (0:floor(dur_sec * op.target_fs)-1)' / op.target_fs;
        y_rs = interp1(t_old, double(y_mono), t_new, 'linear', 0);
    end
else
    y_rs = double(y_mono);
end
clear y_mono;

% Zero-phase bandpass filter within shared speech band
nyq = op.target_fs / 2;
bp = [max(20, op.bandpass_hz(1)), min(0.95 * nyq, op.bandpass_hz(2))];
if exist('butter', 'file') && exist('filtfilt', 'file')
    [b, a] = butter(4, bp / nyq, 'bandpass');
    y_cond = filtfilt(b, a, y_rs);
else
    % FFT brick-wall fallback if Signal Processing Toolbox is unavailable
    N = numel(y_rs);
    Y = fft(y_rs);
    freqs = (0:N-1)' * (op.target_fs / N);
    freqs(freqs > nyq) = op.target_fs - freqs(freqs > nyq);
    mask = (freqs >= bp(1)) & (freqs <= bp(2));
    y_cond = real(ifft(Y .* mask));
end
end

% -------------------------------------------------------------------------
function res = align_audio_pair(y_mas, y_mov, op, landmarks)
fs = op.target_fs;
n_mas = numel(y_mas);
n_mov = numel(y_mov);

% --- STAGE 1: Broadband Speech Envelope Macro-Alignment ---
win_samps = max(1, round(op.env_win_ms / 1000 * fs));
decim_fac = max(1, round(fs / op.env_fs));
fs_env    = fs / decim_fac;

env_mas = sqrt(movmean(y_mas.^2, win_samps));
env_mov = sqrt(movmean(y_mov.^2, win_samps));
env_mas = env_mas(1:decim_fac:end);
env_mov = env_mov(1:decim_fac:end);

% Sliding local z-score to handle speaker distance differences
z_win = max(5, round(op.local_z_sec * fs_env));
env_mas_z = (env_mas - movmean(env_mas, z_win)) ./ max(movstd(env_mas, z_win), 1e-6);
env_mov_z = (env_mov - movmean(env_mov, z_win)) ./ max(movstd(env_mov, z_win), 1e-6);

% Compute FFT normalized cross-correlation across all valid overlaps for QC & unseeded lag
[r_global, lags_env] = norm_xcorr_fft(env_mas_z, env_mov_z);
lags_sec = lags_env / fs_env;
[r_glob_peak, idx_glob] = max(r_global);
lag_glob_sec = parabolic_peak_interp(r_global, lags_sec, idx_glob);

% Peak-to-sidelobe Z-score (exclude +-0.15 s around peak)
bg_mask = abs(lags_sec - lag_glob_sec) > 0.15;
if any(bg_mask)
    mu_bg  = mean(r_global(bg_mask));
    sig_bg = max(std(r_global(bg_mask)), 1e-6);
    peak_z = (r_glob_peak - mu_bg) / sig_bg;
else
    peak_z = 0;
end

if landmarks.valid
    w_seed      = landmarks.w_manual;
    off_seed    = landmarks.off_manual;
    search_rad  = op.fine_search_rad_seeded_sec;
else
    w_seed      = 1.0;
    off_seed    = lag_glob_sec;
    search_rad  = op.fine_search_rad_unseeded_sec;
end

% --- STAGE 2: Piecewise Local Cross-Correlation & Linear Drift Fit ---
fine_smooth = max(1, round(op.fine_env_smooth_ms / 1000 * fs));
sig_mas = movmean(abs(y_mas), fine_smooth);
sig_mov = movmean(abs(y_mov), fine_smooth);

win_n  = round(op.fine_win_sec * fs);
step_n = round(op.fine_step_sec * fs);
max_lag_samps = round(search_rad * fs);

% Candidate window centers across moving recording
half_w = floor(win_n / 2);
centers_mov = (half_w + 1) : step_n : (n_mov - half_w);
if numel(centers_mov) < 6 && n_mov > win_n
    centers_mov = round(linspace(half_w + 1, n_mov - half_w, min(12, floor(n_mov / (2*fs)))));
end

% Energy threshold to skip pure-silence windows
mov_rms_all = rms(y_mov);
active_thresh = 0.15 * mov_rms_all;

t_mov_pts = [];
t_mas_pts = [];
r_win_pts = [];

for k = 1:numel(centers_mov)
    c_mov = centers_mov(k);
    idx1_mov = c_mov - half_w;
    idx2_mov = idx1_mov + win_n - 1;
    seg_mov_raw = y_mov(idx1_mov:idx2_mov);

    if rms(seg_mov_raw) < active_thresh
        continue;
    end

    t_c_mov = (c_mov - 1) / fs;
    t_pred_mas = w_seed * t_c_mov + off_seed;
    c_mas_pred = round(t_pred_mas * fs) + 1;

    idx1_mas = c_mas_pred - half_w - max_lag_samps;
    idx2_mas = idx1_mas + win_n - 1 + 2 * max_lag_samps;
    if idx1_mas < 1 || idx2_mas > n_mas
        continue;
    end

    seg_mov = sig_mov(idx1_mov:idx2_mov);
    seg_mas = sig_mas(idx1_mas:idx2_mas);
    seg_mov = seg_mov - mean(seg_mov);
    seg_mas = seg_mas - mean(seg_mas);

    if std(seg_mov) < 1e-8 || std(seg_mas) < 1e-8
        continue;
    end

    % Valid cross-correlation across [-max_lag_samps, +max_lag_samps]
    n_lags = 2 * max_lag_samps + 1;
    r_loc = zeros(n_lags, 1);
    norm_mov = norm(seg_mov);
    for L = 1:n_lags
        sub_mas = seg_mas(L : L + win_n - 1);
        denom = norm_mov * norm(sub_mas);
        if denom > 0
            r_loc(L) = (seg_mov' * sub_mas) / denom;
        end
    end

    [r_max, idx_max] = max(r_loc);
    lag_grid_sec = (-max_lag_samps : max_lag_samps)' / fs;
    fine_delta_sec = parabolic_peak_interp(r_loc, lag_grid_sec, idx_max);

    t_c_mas = t_pred_mas + fine_delta_sec;
    t_mov_pts(end+1, 1) = t_c_mov; %#ok<AGROW>
    t_mas_pts(end+1, 1) = t_c_mas; %#ok<AGROW>
    r_win_pts(end+1, 1) = r_max;   %#ok<AGROW>
end

% Fit linear clock model: t_mas = w * t_mov + offset
n_active = numel(t_mov_pts);
if n_active >= 2
    [w_fit, off_fit, inlier_mask] = irls_linear_fit(t_mov_pts, t_mas_pts, r_win_pts, op.inlier_tol_ms / 1000);
else
    w_fit = w_seed;
    off_fit = off_seed;
    inlier_mask = false(n_active, 1);
end

n_inliers = sum(inlier_mask);
if n_active > 0
    inlier_frac = n_inliers / n_active;
    r_median = median(r_win_pts);
else
    inlier_frac = 0;
    r_median = 0;
end

if n_inliers >= 2
    resid_sec = t_mas_pts(inlier_mask) - (w_fit * t_mov_pts(inlier_mask) + off_fit);
    resid_rmse_ms = sqrt(mean(resid_sec.^2)) * 1000;
else
    resid_rmse_ms = Inf;
end

clock_drift_ppm = (w_fit - 1) * 1e6;

% Manual landmark discrepancy check
if landmarks.valid
    err_t1_ms = ((w_fit * landmarks.t_mov_1 + off_fit) - landmarks.t_mas_1) * 1000;
    err_t2_ms = ((w_fit * landmarks.t_mov_2 + off_fit) - landmarks.t_mas_2) * 1000;
    max_err_ms = max(abs([err_t1_ms, err_t2_ms]));
else
    err_t1_ms = NaN;
    err_t2_ms = NaN;
    max_err_ms = NaN;
end

% Evaluate audio QC cutoffs
fail_reasons = {};
if r_median < op.qc_min_corr
    fail_reasons{end+1} = sprintf('r_median=%.3f < %.3f', r_median, op.qc_min_corr);
end
if peak_z < op.qc_min_peak_z
    fail_reasons{end+1} = sprintf('peak_z=%.1f < %.1f', peak_z, op.qc_min_peak_z);
end
if inlier_frac < op.qc_min_inlier_frac
    fail_reasons{end+1} = sprintf('inlier_frac=%.2f < %.2f', inlier_frac, op.qc_min_inlier_frac);
end
if n_inliers < op.qc_min_inlier_win
    fail_reasons{end+1} = sprintf('n_inliers=%d < %d', n_inliers, op.qc_min_inlier_win);
end
if resid_rmse_ms > op.qc_max_resid_ms
    fail_reasons{end+1} = sprintf('resid_rmse=%.1fms > %.1fms', resid_rmse_ms, op.qc_max_resid_ms);
end
if abs(clock_drift_ppm) > op.qc_max_warp_ppm
    fail_reasons{end+1} = sprintf('|drift|=%.0fppm > %.0fppm', abs(clock_drift_ppm), op.qc_max_warp_ppm);
end

res = struct( ...
    'warpfactor', w_fit, ...
    'offset', off_fit, ...
    'r_median', r_median, ...
    'peak_z', peak_z, ...
    'inlier_frac', inlier_frac, ...
    'n_inliers', n_inliers, ...
    'n_active', n_active, ...
    'resid_rmse_ms', resid_rmse_ms, ...
    'clock_drift_ppm', clock_drift_ppm, ...
    'landmark_err_t1_ms', err_t1_ms, ...
    'landmark_err_t2_ms', err_t2_ms, ...
    'landmark_max_err_ms', max_err_ms, ...
    'audio_qc_passed', isempty(fail_reasons), ...
    'qc_fail_reasons', {fail_reasons}, ...
    't_mov_pts', t_mov_pts, ...
    't_mas_pts', t_mas_pts, ...
    'r_win_pts', r_win_pts, ...
    'inlier_mask', inlier_mask, ...
    'lags_sec', lags_sec, ...
    'r_global', r_global);
end

% -------------------------------------------------------------------------
function [r_norm, lags] = norm_xcorr_fft(x, y)
% Fast normalized cross-correlation between x (master) and y (moving).
% Positive lag means y occurs later in x (t_x = t_y + lag).
nx = numel(x);
ny = numel(y);
n_fft = 2^nextpow2(nx + ny - 1);

X = fft(x, n_fft);
Y = fft(y, n_fft);
cc = real(ifft(X .* conj(Y)));

% Arrange lags from -(ny-1) to +(nx-1)
cc_ordered = [cc(end - ny + 2 : end); cc(1 : nx)];
lags = (-(ny - 1) : (nx - 1))';

% Exact overlapping window norms via cumulative sums
x2_cs = [0; cumsum(x.^2)];
y2_cs = [0; cumsum(y.^2)];

k = (1:numel(lags))';
x_start = max(1, lags + 1);
x_end   = min(nx, lags + ny);
y_start = max(1, 1 - lags);
y_end   = min(ny, nx - lags);

norm_x = sqrt(max(0, x2_cs(x_end + 1) - x2_cs(x_start)));
norm_y = sqrt(max(0, y2_cs(y_end + 1) - y2_cs(y_start)));
denom = norm_x .* norm_y;

r_norm = zeros(size(cc_ordered));
valid = denom > 1e-8 & (x_end - x_start + 1) >= round(0.1 * min(nx, ny));
r_norm(valid) = cc_ordered(valid) ./ denom(valid);
end

% -------------------------------------------------------------------------
function x_interp = parabolic_peak_interp(y, x, idx)
if idx <= 1 || idx >= numel(y)
    x_interp = x(idx);
    return;
end
y1 = y(idx - 1);
y2 = y(idx);
y3 = y(idx + 1);
denom = (y1 - 2*y2 + y3);
if abs(denom) < 1e-12
    x_interp = x(idx);
else
    delta = 0.5 * (y1 - y3) / denom;
    delta = max(-0.5, min(0.5, delta));
    dx = x(2) - x(1);
    x_interp = x(idx) + delta * dx;
end
end

% -------------------------------------------------------------------------
function [w, off, inliers] = irls_linear_fit(x, y, r_weights, tol_sec)
% Iteratively Reweighted Least Squares fit for y = w * x + off.
n = numel(x);
% Initialize using median slope of pairwise anchor points (Theil-Sen style)
slopes = [];
for i = 1:n-1
    for j = i+1:n
        if abs(x(j) - x(i)) > 5.0
            slopes(end+1, 1) = (y(j) - y(i)) / (x(j) - x(i)); %#ok<AGROW>
        end
    end
end
if isempty(slopes)
    w = 1.0;
else
    w = median(slopes);
end
off = median(y - w * x);

wts = max(r_weights, 0.01).^2;
for iter = 1:25
    resid = y - (w * x + off);
    scale = max(1.4826 * median(abs(resid)), tol_sec / 2);
    u = resid / (4.685 * scale);
    bi_w = zeros(size(u));
    m = abs(u) < 1;
    bi_w(m) = (1 - u(m).^2).^2;
    eff_w = wts .* bi_w;
    if sum(eff_w > 0) < 2
        break;
    end
    W = diag(eff_w);
    A = [x, ones(n, 1)];
    theta = (A' * W * A) \ (A' * W * y);
    w = theta(1);
    off = theta(2);
end

final_resid = abs(y - (w * x + off));
inliers = final_resid <= tol_sec;
if sum(inliers) >= 2
    p = polyfit(x(inliers), y(inliers), 1);
    w = p(1);
    off = p(2);
end
end

% -------------------------------------------------------------------------
function plot_alignment_qc(y_mas, y_mov, res, landmarks, op, g_label, mov_disp, master_disp, fig_path)
vis_str = 'off';
if is_true_flag(op.show_qc_figs)
    vis_str = 'on';
end
fig = figure('Visible', vis_str, 'Position', [100 100 1150 750], 'Color', 'w');

% Panel 1: Piecewise Clock Drift & Linear Fit
subplot(2, 2, 1); hold on;
if ~isempty(res.t_mov_pts)
    lag_ms = (res.t_mas_pts - res.t_mov_pts) * 1000;
    in_m = res.inlier_mask;
    scatter(res.t_mov_pts(in_m)/60, lag_ms(in_m), 36, res.r_win_pts(in_m), 'filled');
    if any(~in_m)
        scatter(res.t_mov_pts(~in_m)/60, lag_ms(~in_m), 36, [0.7 0.7 0.7], 'x', 'LineWidth', 1.2);
    end
    cb = colorbar; ylabel(cb, 'Window r');
    t_line = linspace(min(res.t_mov_pts), max(res.t_mov_pts), 100);
    plot(t_line/60, ((res.warpfactor - 1)*t_line + res.offset)*1000, 'r-', 'LineWidth', 1.6);
end
if landmarks.valid
    plot([landmarks.t_mov_1, landmarks.t_mov_2]/60, ...
         [(landmarks.t_mas_1 - landmarks.t_mov_1)*1000, (landmarks.t_mas_2 - landmarks.t_mov_2)*1000], ...
         'ko--', 'MarkerSize', 8, 'LineWidth', 1.3);
end
xlabel('Moving Time (min)');
ylabel('Lag t_{mas} - t_{mov} (ms)');
title(sprintf('Clock Drift Fit: w=%.7f (%+.1f ppm), off=%.3fs', ...
    res.warpfactor, res.clock_drift_ppm, res.offset));
grid on;

% Panel 2: Stage 1 Global Envelope Cross-Correlation
subplot(2, 2, 2);
plot(res.lags_sec, res.r_global, 'b-'); hold on;
xline(res.offset, 'r--', 'LineWidth', 1.3);
xlabel('Lag (s)');
ylabel('Normalized Envelope r');
title(sprintf('Macro Envelope XCorr (Z_{peak}=%.1f, r_{med}=%.2f)', res.peak_z, res.r_median));
grid on;

% Panel 3 & 4: Zoomed 2s Envelope Overlays at Start and End of Overlap
fs = op.target_fs;
env_win = round(0.015 * fs);
e_mas = movmean(abs(y_mas), env_win);
e_mov = movmean(abs(y_mov), env_win);
t_mas = (0:numel(e_mas)-1)' / fs;
t_mov_warped = res.warpfactor * ((0:numel(e_mov)-1)' / fs) + res.offset;

ov_start = max(t_mas(1), t_mov_warped(1)) + 1.0;
ov_end   = min(t_mas(end), t_mov_warped(end)) - 3.0;
if ov_end <= ov_start
    ov_start = max(t_mas(1), t_mov_warped(1));
    ov_end   = ov_start + 2.0;
end
zoom_starts = [ov_start + 0.15*(ov_end - ov_start), ov_start + 0.80*(ov_end - ov_start)];
zoom_titles = {'Aligned Envelope Zoom (Early Overlap)', 'Aligned Envelope Zoom (Late Overlap)'};

for z = 1:2
    subplot(2, 2, 2 + z); hold on;
    z1 = zoom_starts(z); z2 = z1 + 2.0;
    m_idx = t_mas >= z1 & t_mas <= z2;
    v_idx = t_mov_warped >= z1 & t_mov_warped <= z2;
    if any(m_idx)
        em = e_mas(m_idx) / max(max(e_mas(m_idx)), 1e-6);
        plot(t_mas(m_idx), em, 'k-', 'LineWidth', 1.1);
    end
    if any(v_idx)
        ev = e_mov(v_idx) / max(max(e_mov(v_idx)), 1e-6);
        plot(t_mov_warped(v_idx), -ev, 'b-', 'LineWidth', 1.1);
    end
    yline(0, 'k:');
    xlabel('Master Time (s)');
    ylabel('Norm Env (+Mas / -Mov)');
    title(zoom_titles{z});
    grid on;
end

sgtitle(sprintf('[%s] %s -> %s', g_label, mov_disp, master_disp), 'Interpreter', 'none');

if op.qc_save_figs
    exportgraphics(fig, fig_path, 'Resolution', 150);
end
if ~is_true_flag(op.show_qc_figs)
    close(fig);
end
end

% -------------------------------------------------------------------------
function [synced_annot_paths, total_clips] = process_annotation_files(mov_audio_path, orig_media_path, warpfactor, offset, op)
synced_annot_paths = {};
total_clips = 0;

[mov_dir, mov_stem, ~] = fileparts(mov_audio_path);
[~, orig_stem, ~] = fileparts(orig_media_path);
stems = unique({mov_stem, orig_stem}, 'stable');

% Determine search directories
search_dirs = op.annot_dir;
if isempty(search_dirs)
    search_dirs = {mov_dir};
end

% Build candidate annotation filename stems
cand_stems = {};
for s = 1:numel(stems)
    st = stems{s};
    for k = 1:numel(op.annotation_file_suffixes)
        suf = strtrim(char(op.annotation_file_suffixes{k}));
        if isempty(suf)
            cand_stems{end+1} = st; %#ok<AGROW>
        else
            cand_stems{end+1} = [st, suf];       %#ok<AGROW>
            cand_stems{end+1} = [st, '_', suf];  %#ok<AGROW>
            cand_stems{end+1} = [st, '-', suf];  %#ok<AGROW>
        end
    end
end
cand_stems = unique(cand_stems, 'stable');
valid_exts = {'.xlsx', '.tsv'};

matched_files = {};
matched_base_dirs = {};
matched_in_unsynced = [];

for d = 1:numel(search_dirs)
    base_d = char(search_dirs{d});
    if ~isfolder(base_d), continue; end
    unsync_d = fullfile(base_d, 'unsynced');

    sub_checks = {base_d, false; unsync_d, true};
    for sc = 1:size(sub_checks, 1)
        chk_dir = sub_checks{sc, 1};
        is_unsync = sub_checks{sc, 2};
        if ~isfolder(chk_dir), continue; end

        listing = dir(chk_dir);
        for f = 1:numel(listing)
            if listing(f).isdir, continue; end
            [~, f_stem, f_ext] = fileparts(listing(f).name);
            if ~any(strcmpi(f_ext, valid_exts)), continue; end
            % Skip already-synced files
            if endsWith(f_stem, '_sync', 'IgnoreCase', true), continue; end

            if any(strcmpi(f_stem, cand_stems))
                full_f = fullfile(chk_dir, listing(f).name);
                if ~any(strcmp(matched_files, full_f))
                    matched_files{end+1} = full_f; %#ok<AGROW>
                    matched_base_dirs{end+1} = base_d; %#ok<AGROW>
                    matched_in_unsynced(end+1) = is_unsync; %#ok<AGROW>
                end
            end
        end
    end
end

for m = 1:numel(matched_files)
    orig_f = matched_files{m};
    base_d = matched_base_dirs{m};
    was_in_unsynced = matched_in_unsynced(m);
    [~, orig_f_stem, orig_f_ext] = fileparts(orig_f);

    % Read annotation table
    if strcmpi(orig_f_ext, '.tsv')
        atbl = readtable(orig_f, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
    else
        atbl = readtable(orig_f, 'TextType', 'string');
    end

    if ~all(ismember({'starts_file_relative', 'ends_file_relative'}, atbl.Properties.VariableNames))
        warning('sync_audio_and_annot:MissingAnnotCols', ...
            'Annotation file %s lacks starts_file_relative/ends_file_relative; skipping.', orig_f);
        continue;
    end

    s_rel = get_numeric_col(atbl, 'starts_file_relative');
    e_rel = get_numeric_col(atbl, 'ends_file_relative');

    starts_mas = warpfactor * s_rel + offset;
    ends_mas   = warpfactor * e_rel + offset;
    dur_mas    = ends_mas - starts_mas;

    % Remove old starts, ends, duration, starts_file_relative, ends_file_relative
    drop_cols = intersect(atbl.Properties.VariableNames, ...
        {'starts', 'ends', 'duration', 'starts_file_relative', 'ends_file_relative'}, 'stable');
    other_cols = setdiff(atbl.Properties.VariableNames, drop_cols, 'stable');

    new_atbl = table(starts_mas, ends_mas, dur_mas, ...
        'VariableNames', {'starts', 'ends', 'duration'});
    if ~isempty(other_cols)
        new_atbl = [new_atbl, atbl(:, other_cols)]; %#ok<AGROW>
    end
    new_atbl.starts_file_relative = s_rel;
    new_atbl.ends_file_relative   = e_rel;

    % Move original file into base_d/unsynced if not already there
    unsync_d = fullfile(base_d, 'unsynced');
    if ~isfolder(unsync_d)
        mkdir(unsync_d);
    end
    if ~was_in_unsynced
        dest_unsync = fullfile(unsync_d, [orig_f_stem, orig_f_ext]);
        movefile(orig_f, dest_unsync);
    end

    % Save synced annotation file one level above 'unsynced' (in base_d)
    synced_f = fullfile(base_d, [orig_f_stem, '_sync', orig_f_ext]);
    if strcmpi(orig_f_ext, '.tsv')
        writetable(new_atbl, synced_f, 'FileType', 'text', 'Delimiter', '\t');
    else
        writetable(new_atbl, synced_f);
    end
    synced_annot_paths{end+1} = synced_f; %#ok<AGROW>

    % Chop event-wise audio clips if eventwise_audvid_dir is specified
    if ~isempty(op.eventwise_audvid_dir)
        n_c = chop_eventwise_audio(mov_audio_path, s_rel, e_rel, orig_f_stem, op);
        total_clips = total_clips + n_c;
    end
end
end

% -------------------------------------------------------------------------
function n_written = chop_eventwise_audio(mov_audio_path, s_rel, e_rel, annot_stem, op)
n_written = 0;
ainfo = audioinfo(mov_audio_path);
fs = ainfo.SampleRate;
total_samps = ainfo.TotalSamples;
buf = op.event_audvid_buffer_sec;

out_dir = fullfile(op.eventwise_audvid_dir, annot_stem);
if ~isfolder(out_dir)
    mkdir(out_dir);
end

for r = 1:numel(s_rel)
    if ~isfinite(s_rel(r)) || ~isfinite(e_rel(r))
        continue;
    end
    t1 = max(0, s_rel(r) - buf(1));
    t2 = min(ainfo.Duration, e_rel(r) + buf(2));
    if t2 <= t1
        continue;
    end
    s1 = max(1, floor(t1 * fs) + 1);
    s2 = min(total_samps, ceil(t2 * fs));
    if s2 <= s1
        continue;
    end
    [y_clip, ~] = audioread(mov_audio_path, [s1, s2]);
    clip_name = sprintf('%s_event-%03d.wav', annot_stem, r);
    audiowrite(fullfile(out_dir, clip_name), y_clip, fs);
    n_written = n_written + 1;
end
end

% -------------------------------------------------------------------------
function print_sync_summary(aligned_log, skipped_log, op)
fprintf('\n====================================================================\n');
fprintf('SYNC SUMMARY\n');
fprintf('====================================================================\n');

fprintf('\n1. ALIGNED MOVING -> MASTER PAIRS (%d succeeded):\n', numel(aligned_log));
if isempty(aligned_log)
    fprintf('   (None)\n');
else
    for i = 1:numel(aligned_log)
        a = aligned_log(i);
        fprintf('   [%s] %s  -->  %s\n', a.group, a.moving_file, a.master_file);
        fprintf('      Status: SUCCESS (%s) | warpfactor = %.8f | offset = %.4f s\n', ...
            a.method, a.warpfactor, a.offset);
        fprintf('      QC Metrics: r_median = %.3f | Clock Drift = %+.1f ppm | Landmark Err: dt1 = %+.1f ms, dt2 = %+.1f ms\n', ...
            a.r_median, a.clock_drift_ppm, a.err_t1_ms, a.err_t2_ms);
        if isempty(a.annot_files)
            fprintf('      Synced Annotation Files: (None found)\n');
        else
            for j = 1:numel(a.annot_files)
                fprintf('      Synced Annotation File : %s\n', a.annot_files{j});
            end
            if ~isempty(op.eventwise_audvid_dir)
                fprintf('      Event-wise Audio Clips : %d clip(s) saved under %s\n', ...
                    a.n_clips, op.eventwise_audvid_dir);
            end
        end
    end
end

fprintf('\n2. SKIPPED / FAILED FILES (%d):\n', numel(skipped_log));
if isempty(skipped_log)
    fprintf('   (None)\n');
else
    for i = 1:numel(skipped_log)
        s = skipped_log(i);
        fprintf('   [%s] %s  -->  %s\n', s.group, s.file, s.reason);
    end
end

fprintf('\n3. OUTPUT LOCATIONS:\n');
fprintf('   Synced table saved to : %s\n', op.sync_table_savename);
if op.qc_save_figs
    fprintf('   QC figures saved to   : %s\n', op.qc_figs_savedir);
end
if ~isempty(op.eventwise_audvid_dir)
    fprintf('   Event-wise audio dir  : %s\n', op.eventwise_audvid_dir);
end
fprintf('====================================================================\n\n');
end

% -------------------------------------------------------------------------
function v = get_numeric_col(tbl, colname)
raw = tbl.(colname);
if isnumeric(raw) || islogical(raw)
    v = double(raw);
elseif isstring(raw) || iscellstr(raw) || iscategorical(raw)
    v = str2double(string(raw));
elseif iscell(raw)
    v = nan(numel(raw), 1);
    for i = 1:numel(raw)
        if isnumeric(raw{i}) && ~isempty(raw{i})
            v(i) = double(raw{i});
        elseif ischar(raw{i}) || isstring(raw{i})
            v(i) = str2double(string(raw{i}));
        end
    end
else
    v = nan(height(tbl), 1);
end
end

% -------------------------------------------------------------------------
function c = get_str_col(tbl, colname)
raw = tbl.(colname);
if iscellstr(raw)
    c = raw;
elseif isstring(raw) || iscategorical(raw)
    c = cellstr(string(raw));
elseif isnumeric(raw) || islogical(raw)
    c = cellstr(string(raw));
elseif iscell(raw)
    c = cell(numel(raw), 1);
    for i = 1:numel(raw)
        if isempty(raw{i})
            c{i} = '';
        else
            c{i} = char(string(raw{i}));
        end
    end
else
    c = repmat({''}, height(tbl), 1);
end
c(cellfun(@(x) strcmpi(x, 'missing') || strcmpi(x, '<missing>'), c)) = {''};
end

% -------------------------------------------------------------------------
function s = get_cell_str_at(tbl, colname, idx)
c = get_str_col(tbl, colname);
s = c{idx};
end

% -------------------------------------------------------------------------
function disp_name = get_display_filename(tbl, idx)
vnames = tbl.Properties.VariableNames;
if ismember('filename', vnames)
    disp_name = get_cell_str_at(tbl, 'filename', idx);
elseif ismember('filepath', vnames)
    [~, n, ext] = fileparts(get_cell_str_at(tbl, 'filepath', idx));
    disp_name = [n, ext];
else
    disp_name = sprintf('row_%d', idx);
end
end

% -------------------------------------------------------------------------
function tf = is_true_flag(val)
if islogical(val) || isnumeric(val)
    tf = logical(val);
elseif ischar(val) || isstring(val)
    tf = any(strcmpi(strtrim(char(val)), {'on', 'true', '1', 'yes'}));
else
    tf = false;
end
end