function annotate_manual_ephys_artifacts(cfg)
%ANNOTATE_MANUAL_EPHYS_ARTIFACTS  Interactive GUI for manual annotation of
%   intra-operative ephys artifacts in a FieldTrip raw structure.
%
% Usage
%   annotate_manual_ephys_artifacts()
%   annotate_manual_ephys_artifacts(cfg)
%
% cfg fields (all optional)
%   .sub                         – Subject ID string, e.g. 'DM1005'  (prompted if absent)
%   .task                        – Task string name, e.g. 'smsl'                ['smsl']
%   .n_chans_raster              – Channels shown per screen in raster view     [40]
%   .n_chans_timecourse          – Channels shown per screen in timecourse view [20]
%   .iqr_thr                     – Interquartile range threshold for artifacts  [3]
%   .max_timecourse_sample_rate  - Target max downsample display rate (Hz)      [200]
%
% Input file
%   Y:\DBS\derivatives\sub-<SUB>\fieldtrip\sub-<SUB>_ses-intraop_task-<TASK>_ft-raw.mat
%   Must contain a FieldTrip raw struct with fields: trial, time, label.
%
% Output .tsv columns
%   id       – integer row index
%   starts   – window start  [global time / s]
%   ends     – window end    [global time / s]
%   duration – ends - starts [s]
%   label    – channel label matching D.label
%
% Requirements: MATLAB R2018b+ (yline, string arrays, datetime arithmetic)
%==========================================================================
%% 0.  Defaults
%==========================================================================
if nargin < 1 || isempty(cfg), cfg = struct(); end
if ~isfield(cfg,'task'),                        cfg.task                        = 'smsl'; end
if ~isfield(cfg,'n_chans_raster'),              cfg.n_chans_raster              = 40; end
if ~isfield(cfg,'n_chans_timecourse'),          cfg.n_chans_timecourse          = 20; end
if ~isfield(cfg,'iqr_thr'),                     cfg.iqr_thr                     = 3; end
if ~isfield(cfg,'max_timecourse_sample_rate'),  cfg.max_timecourse_sample_rate  = 200; end

%==========================================================================
%% 1.  Subject ID & Task Dialog
%==========================================================================
if ~isfield(cfg,'sub') || isempty(cfg.sub) || ~isfield(cfg,'task') || isempty(cfg.task)
    prompt = {'Enter subject ID:', 'Enter task:', 'Max Timecourse Sample Rate (Hz):'};
    dlgtitle = 'Setup Configuration';
    dims = [1 35];
    definput = {'DM10', 'smsl', num2str(cfg.max_timecourse_sample_rate)};
    ans_ = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(ans_), return; end
    cfg.sub  = strtrim(ans_{1});
    cfg.task = strtrim(ans_{2});
    cfg.max_timecourse_sample_rate = str2double(strtrim(ans_{3}));
end
%==========================================================================
%% 2.  Load FieldTrip file
%==========================================================================
ft_dir = fullfile('Y:\DBS', 'derivatives', ['sub-' cfg.sub], 'fieldtrip');
ft_name = sprintf('sub-%s_ses-intraop_task-%s_ft-raw.mat', cfg.sub, cfg.task);
ft_path = fullfile(ft_dir, ft_name);

fprintf('\n=== FieldTrip file ===\n');
fprintf('  Target Path : %s\n', ft_path);

% Fallback selector if target file path does not exist
if ~exist(ft_path, 'file')
    filter_pattern = sprintf('*task-%s_ft-raw.mat', cfg.task);
    default_search = fullfile(ft_dir, filter_pattern);
    
    [f, p] = uigetfile({'*.mat', 'FieldTrip Data (*.mat)'}, ...
        'Failed to find fieldtrip file; select fieldtrip file', default_search);
    if isequal(f, 0)
        fprintf('Loading cancelled by user.\n');
        return; 
    end
    ft_path = fullfile(p, f);
end

dinfo = dir(ft_path);
if isempty(dinfo)
    errordlg(sprintf('File not found:\n%s', ft_path), 'File Not Found');
    return;
end
fprintf('  Loading Path: %s\n', ft_path);
fprintf('  Size        : %.2f MB\n', dinfo.bytes/1e6);

t_load0 = datetime('now');
fprintf('  Load started: %s\n', char(t_load0,'HH:mm:ss.SSS'));
loaded = load(ft_path);
t_load1 = datetime('now');
fprintf('  Load done   : %s  (%.2f s)\n\n', ...
    char(t_load1,'HH:mm:ss.SSS'), seconds(t_load1 - t_load0));

% Identify FieldTrip raw struct
D = [];
for fv_ = fieldnames(loaded)'
    cand_ = loaded.(fv_{1});
    if isstruct(cand_) && all(isfield(cand_,{'trial','label','time'}))
        D = cand_; break;
    end
end
if isempty(D)
    errordlg('No FieldTrip raw struct (trial/label/time) found in file.', 'Load Error');
    return;
end
data_mat  = D.trial{1};         % [nCh x nTime]
time_vec  = D.time{1};          % [1   x nTime]
ch_labels = cellstr(D.label);   % {nCh x 1}
n_chans   = size(data_mat,1);

%==========================================================================
%% 3.  Shared state  (accessed/modified by nested functions)
%==========================================================================
op = struct();
op.max_timecourse_sample_rate = cfg.max_timecourse_sample_rate;

viewmode  = 'timecourse';    % 'raster' | 'timecourse'
cur_chunk = 1;
% Custom zoom limits
custom_xlim = [];
custom_ylim = [];
% Trace Scaling values
raw_scale_val = 1.0;
clean_scale_val = 1.0;
notch_scale_val = 1.0;

% Processed data targets
data_mat_clean = [];
data_mat_notch = [];

% Artifact table
artifact = table( ...
    zeros(0,1,'double'), zeros(0,1,'double'), ...
    zeros(0,1,'double'), zeros(0,1,'double'), strings(0,1), ...
    'VariableNames',{'id','starts','ends','duration','label'});
% Current gray selections – struct array
SEL = struct('t_start',{},'t_end',{},'ch_start',{},'ch_end',{});
% Rubber-band drag state
is_dragging = false;
drag_x0     = NaN;
drag_y0     = NaN;
rband_h     = [];   % graphics handle, [] = none
% Colormaps (15)
CMAPS = {'parula','turbo','jet','hsv','hot','cool','spring','summer', ...
         'autumn','winter','gray','bone','copper','pink','colorcube'};
cur_cmap   = 'parula';

%==========================================================================
%% 4.  Pre-compute Initial Data & Display Variables
%==========================================================================
fprintf('  Running initial data processing...\n');
cleaning_func_cfg_out = compute_processed_traces(true);

%==========================================================================
%% 5.  Build figure & controls
%==========================================================================
LP_W = 0.19;   % fractional figure width for left panel
fig = figure( ...
    'Name',          ['Artifact Annotator  |  sub-' cfg.sub '  |  task-' cfg.task], ...
    'NumberTitle',   'off', ...
    'Units',         'normalized', ...
    'OuterPosition', [0.01 0.03 0.98 0.94], ...
    'Color',         [0.93 0.93 0.93], ...
    'WindowButtonDownFcn',   @btn_down_cb, ...
    'WindowButtonMotionFcn', @btn_motion_cb, ...
    'WindowButtonUpFcn',     @btn_up_cb, ...
    'WindowKeyPressFcn',     @key_press_cb, ...
    'WindowScrollWheelFcn',  @scroll_zoom_cb); % Interactive mouse wheel integration

% ── Left panel ────────────────────────────────────────────────────────────
lp = uipanel(fig, ...
    'Title','Controls','Units','normalized', ...
    'Position',[0 0 LP_W 1], ...
    'BackgroundColor',[0.87 0.87 0.87],'FontSize',9);
% Stacking helpers
ny_ = 0.985;  ch_ = 0.028;  dh_ = 0.032;

% View mode Toggle (Single Button using HTML for diverse fonts)
yp_ = ny_-ch_; ny_ = ny_-dh_;
btn_view = uicontrol(lp,'Style','pushbutton', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'Callback',@(~,~) viewmode_cb());

% Fixed Legend Control (Split into three text boxes to avoid HTML completely)
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','Black=Raw', ...
    'Units','normalized','Position',[0.02 yp_ 0.30 ch_], ...
    'HorizontalAlignment','center','BackgroundColor',[0.87 0.87 0.87], ...
    'FontWeight','bold', 'ForegroundColor','k');
uicontrol(lp,'Style','text','String','Blue=Clean', ...
    'Units','normalized','Position',[0.33 yp_ 0.30 ch_], ...
    'HorizontalAlignment','center','BackgroundColor',[0.87 0.87 0.87], ...
    'FontWeight','bold', 'ForegroundColor',[0 0 0.8]);
uicontrol(lp,'Style','text','String','Red=Notched', ...
    'Units','normalized','Position',[0.64 yp_ 0.34 ch_], ...
    'HorizontalAlignment','center','BackgroundColor',[0.87 0.87 0.87], ...
    'FontWeight','bold', 'ForegroundColor',[0.6 0 0]);

% Display Sample Rate Sub-Legend Text Box
yp_ = ny_-ch_; ny_ = ny_-dh_;
txt_fs_display = uicontrol(lp,'Style','text','String','Sample Rate: -- Hz', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87], ...
    'FontWeight','bold', 'ForegroundColor','k');

% Colormap Selector
yp_ = ny_-ch_; ny_ = ny_-dh_;
lbl_cmap = uicontrol(lp,'Style','text','String','Colormap:', ...
    'Units','normalized','Position',[0.05 yp_ 0.40 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
dd_cmap = uicontrol(lp,'Style','popupmenu','String',CMAPS, ...
    'Units','normalized','Position',[0.45 yp_ 0.50 ch_],'Callback',@cmap_cb);

% Channel Visibility Configurations
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','N Visible Chans Raster:', ...
    'Units','normalized','Position',[0.05 yp_ 0.50 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
ed_raster = uicontrol(lp,'Style','edit','String',num2str(cfg.n_chans_raster), ...
    'Units','normalized','Position',[0.60 yp_ 0.35 ch_]);
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','N Visible Chans Timecourse:', ...
    'Units','normalized','Position',[0.05 yp_ 0.50 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
ed_timecourse = uicontrol(lp,'Style','edit','String',num2str(cfg.n_chans_timecourse), ...
    'Units','normalized','Position',[0.60 yp_ 0.35 ch_]);
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Update Chans', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@update_chans_cb);
ny_ = ny_ - 0.010;   % spacer

% Trace Scaling & Visibility Controls
% RAW
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','Raw scale:', ...
    'Units','normalized','Position',[0.02 yp_ 0.38 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
ed_raw_scale = uicontrol(lp,'Style','text','String',num2str(raw_scale_val,'%.2f'), ...
    'Units','normalized','Position',[0.40 yp_ 0.15 ch_], 'BackgroundColor','w');
uicontrol(lp,'Style','pushbutton','String','▲', ...
    'Units','normalized','Position',[0.57 yp_ 0.12 ch_],'Callback',@(~,~) scale_cb('raw', 0.25));
uicontrol(lp,'Style','pushbutton','String','▼', ...
    'Units','normalized','Position',[0.70 yp_ 0.12 ch_],'Callback',@(~,~) scale_cb('raw', -0.25));
btn_vis_raw = uicontrol(lp,'Style','pushbutton','String','✓', ...
    'Units','normalized','Position',[0.84 yp_ 0.12 ch_],'Callback',@toggle_vis_cb);

% CLEANED
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','Cleaned scale:', ...
    'Units','normalized','Position',[0.02 yp_ 0.38 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
ed_clean_scale = uicontrol(lp,'Style','text','String',num2str(clean_scale_val,'%.2f'), ...
    'Units','normalized','Position',[0.40 yp_ 0.15 ch_], 'BackgroundColor','w');
uicontrol(lp,'Style','pushbutton','String','▲', ...
    'Units','normalized','Position',[0.57 yp_ 0.12 ch_],'Callback',@(~,~) scale_cb('clean', 0.25));
uicontrol(lp,'Style','pushbutton','String','▼', ...
    'Units','normalized','Position',[0.70 yp_ 0.12 ch_],'Callback',@(~,~) scale_cb('clean', -0.25));
btn_vis_clean = uicontrol(lp,'Style','pushbutton','String','✓', ...
    'Units','normalized','Position',[0.84 yp_ 0.12 ch_],'Callback',@toggle_vis_cb);

% NOTCHED
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','Notch scale:', ...
    'Units','normalized','Position',[0.02 yp_ 0.38 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
ed_notch_scale = uicontrol(lp,'Style','text','String',num2str(notch_scale_val,'%.2f'), ...
    'Units','normalized','Position',[0.40 yp_ 0.15 ch_], 'BackgroundColor','w');
uicontrol(lp,'Style','pushbutton','String','▲', ...
    'Units','normalized','Position',[0.57 yp_ 0.12 ch_],'Callback',@(~,~) scale_cb('notch', 0.25));
uicontrol(lp,'Style','pushbutton','String','▼', ...
    'Units','normalized','Position',[0.70 yp_ 0.12 ch_],'Callback',@(~,~) scale_cb('notch', -0.25));
btn_vis_notch = uicontrol(lp,'Style','pushbutton','String','✓', ...
    'Units','normalized','Position',[0.84 yp_ 0.12 ch_],'Callback',@toggle_vis_cb);
ny_ = ny_ - 0.010;   % spacer

% Channel groups navigation
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','Channel group:', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
yp_ = ny_-ch_; ny_ = ny_-dh_;
dd_chunk = uicontrol(lp,'Style','popupmenu','String',{'...'}, ...
    'Units','normalized','Position',[0.05 yp_ 0.65 ch_],'Callback',@chunk_cb);
uicontrol(lp,'Style','pushbutton','String','▲', ...
    'Units','normalized','Position',[0.72 yp_ 0.10 ch_],'Callback',@chunk_prev_cb);
uicontrol(lp,'Style','pushbutton','String','▼', ...
    'Units','normalized','Position',[0.84 yp_ 0.10 ch_],'Callback',@chunk_next_cb);
ny_ = ny_ - 0.018;   % spacer

% Zoom controls
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Zoom selection (q)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'Callback',@zoom_sel_cb);
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Zoom full (w)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'Callback',@zoom_full_cb);
ny_ = ny_ - 0.018;   % spacer

% Artifact actions
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Add artifact (a)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'FontWeight','bold','Callback',@add_artifact_cb);
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Remove selected artifact (r)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@remove_artifact_cb);
ny_ = ny_ - 0.018;   % spacer

% File Saving / Loading
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Load artifact mask (o)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@load_artifact_cb);
yp_ = ny_-ch_; 
uicontrol(lp,'Style','pushbutton','String','Save artifact mask (k)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@save_artifact_cb);

% Sub-panel for Cleaning Configuration
pnl_clean = uipanel(lp, 'Title','Cleaning Config', 'Units','normalized', ...
    'Position',[0.02 0.01 0.96 0.23], 'BackgroundColor',[0.87 0.87 0.87]);

tt_spike = 'estimated duration of a spike artifact in seconds, from spike onset to peak';
tt_iqr   = 'threshold to identify outliers (e.g. outlier > 75th percentile + iqr_thr*interquartile range)';
tt_fc    = 'Cutoff frequency for high-pass filter';

uicontrol(pnl_clean,'Style','text','String','Expected spike dur (s):', ...
    'Units','normalized','Position',[0.02 0.76 0.65 0.20], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87], ...
    'TooltipString', tt_spike);
ed_spike_dur = uicontrol(pnl_clean,'Style','edit', ...
    'String',num2str(cleaning_func_cfg_out.spike_dur), ...
    'Units','normalized','Position',[0.70 0.76 0.25 0.20], ...
    'TooltipString', tt_spike);

uicontrol(pnl_clean,'Style','text','String','Outlier IQR threshold:', ...
    'Units','normalized','Position',[0.02 0.52 0.65 0.20], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87], ...
    'TooltipString', tt_iqr);
ed_iqr_thr = uicontrol(pnl_clean,'Style','edit', ...
    'String',num2str(cleaning_func_cfg_out.iqr_thr), ...
    'Units','normalized','Position',[0.70 0.52 0.25 0.20], ...
    'TooltipString', tt_iqr);

init_fc = '0';
if isfield(cleaning_func_cfg_out, 'f_c')
    init_fc = num2str(cleaning_func_cfg_out.f_c);
end

uicontrol(pnl_clean,'Style','text','String','High pass cutoff (Hz):', ...
    'Units','normalized','Position',[0.02 0.28 0.65 0.20], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87], ...
    'TooltipString', tt_fc);
ed_fc = uicontrol(pnl_clean,'Style','edit', ...
    'String',init_fc, ...
    'Units','normalized','Position',[0.70 0.28 0.25 0.20], ...
    'TooltipString', tt_fc);

uicontrol(pnl_clean,'Style','pushbutton','String','Update traces', ...
    'Units','normalized','Position',[0.05 0.04 0.45 0.18], ...
    'Callback',@(~,~) compute_processed_traces(false));

chk_manual_mask = uicontrol(pnl_clean,'Style','checkbox','String','Include masks', ...
    'Units','normalized','Position',[0.52 0.04 0.45 0.18], ...
    'BackgroundColor',[0.87 0.87 0.87], 'Value', 0);

% Data Selection Info Banner
txt_sel_info = uicontrol(fig, 'Style','text', 'Units','normalized', ...
    'Position',[LP_W+0.04, 0.02, 0.6, 0.04], 'HorizontalAlignment','left', ...
    'BackgroundColor', [0.93 0.93 0.93], 'FontSize', 10, ...
    'String', 'Selection: None');

% ── Main axes ─────────────────────────────────────────────────────────────
ax = axes(fig,'Units','normalized', ...
    'Position',[LP_W+0.04, 0.09, 1-LP_W-0.06, 0.87]);

%==========================================================================
%% 6.  Initialise display
%==========================================================================
% Auto-load existing artifacts if present
annot_dir = fullfile('Y:\DBS','derivatives',['sub-' cfg.sub],'annot');
auto_art_file = sprintf('sub-%s_ses-intraop_task-%s_artifact-manual.tsv', cfg.sub, cfg.task);
auto_art_path = fullfile(annot_dir, auto_art_file);
if exist(auto_art_path, 'file')
    fprintf('  Auto-loading manual artifact file...\n');
    load_artifact_file(auto_art_path);
end

update_viewmode_button();
update_chunk_dropdown();
update_plot();

%==========================================================================
%%   ═══════════  N E S T E D   F U N C T I O N S  ═══════════
%==========================================================================

% ── Process Filtering (Clean + Notch) ─────────────────────────────────────
    function [cfg_out] = compute_processed_traces(is_init)
        if nargin < 1, is_init = false; end
        
        if ~is_init
            set(fig, 'pointer', 'watch'); drawnow;
        end
        
        dlg = dialog('Name', 'Please wait', 'Position', [500 500 250 80]);
        txt_dlg = uicontrol(dlg, 'Style', 'text', 'Position', [20 20 210 40], ...
            'String', 'Updating traces... 0.0 s', 'FontSize', 10);
        t0 = tic;
        tmr = timer('ExecutionMode','fixedSpacing', 'Period',0.1, ...
                    'TimerFcn', @(~,~) update_stopwatch_dlg(txt_dlg, t0));
        start(tmr);

        try
            cfg_clean = struct();
            if ~is_init
                cfg_clean.spike_dur = str2double(ed_spike_dur.String);
                cfg_clean.iqr_thr = str2double(ed_iqr_thr.String);
                cfg_clean.f_c = str2double(ed_fc.String);
                if chk_manual_mask.Value == 1
                    cfg_clean.add_mask = true;
                    cfg_clean.mask_table = artifact;
                else
                    cfg_clean.add_mask = false;
                end
            end
            
            % 1. Cleaned & HPF
            [D_cleaned_new, cfg_out] = hpf_and_instantaneous_artifact_mask(D, cfg_clean);
            data_mat_clean = D_cleaned_new.trial{1};
            % Safeguard against dimension transposition bugs
            if size(data_mat_clean, 2) < size(data_mat_clean, 1)
                data_mat_clean = data_mat_clean';
            end
            
            % 2. Notch Harmonics Filtering
            cfg_notch_filt = [];
            [D_notched_new, ~] = notch_harmonics_filter(D_cleaned_new, cfg_notch_filt);
            data_mat_notch = D_notched_new.trial{1};
            % Safeguard against dimension transposition bugs
            if size(data_mat_notch, 2) < size(data_mat_notch, 1)
                data_mat_notch = data_mat_notch';
            end
            
        catch ME
            stop(tmr); delete(tmr);
            if isgraphics(dlg), close(dlg); end
            if ~is_init, set(fig, 'pointer', 'arrow'); end
            errordlg(sprintf('Processing failed:\n%s', ME.message), 'Processing Error');
            cfg_out = [];
            return;
        end
        
        stop(tmr); delete(tmr);
        if isgraphics(dlg), close(dlg); end
        
        if ~is_init
            set(fig, 'pointer', 'arrow');
            ed_spike_dur.String = num2str(cfg_out.spike_dur);
            ed_iqr_thr.String = num2str(cfg_out.iqr_thr);
            if isfield(cfg_out, 'f_c')
                ed_fc.String = num2str(cfg_out.f_c);
            end
            update_plot();
        end
    end

    function update_stopwatch_dlg(txt_dlg, t0)
        if isgraphics(txt_dlg)
            el = toc(t0);
            txt_dlg.String = sprintf('Updating traces...\nElapsed time: %.1f s', el);
        end
    end

% ── Mouse Scroll Wheel Zoom Callback ──────────────────────────────────────
    function scroll_zoom_cb(~, event)
        cp = ax_coords_clamped();
        if isempty(cp), return; end 
        
        if isempty(custom_xlim)
            xl = [time_vec(1), time_vec(end)];
        else
            xl = custom_xlim;
        end
        
        width = diff(xl);
        mouse_x = cp(1);
        pct = (mouse_x - xl(1)) / width; 
        
        zoom_factor = 0.15; 
        if event.VerticalScrollCount > 0
            new_width = width * (1 + zoom_factor); 
        else
            new_width = width * (1 - zoom_factor); 
        end
        
        total_dur = time_vec(end) - time_vec(1);
        if new_width > total_dur * 5, new_width = total_dur * 5; end
        if new_width < 0.01, new_width = 0.01; end
        
        custom_xlim = [mouse_x - pct * new_width, mouse_x + (1 - pct) * new_width];
        update_plot();
    end

% ── Simple Navigation Helpers ─────────────────────────────────────────────
    function n = n_per_chunk()
        if strcmp(viewmode,'raster'), n = cfg.n_chans_raster;
        else,                         n = cfg.n_chans_timecourse; end
    end

    function nc = n_chunks()
        nc = ceil(n_chans / n_per_chunk());
    end

    function vis = get_vis_chans()
        npc = n_per_chunk();
        s   = (cur_chunk-1)*npc + 1;
        e   = min(s+npc-1, n_chans);
        vis = s:e;
    end

    function gi = lbl2idx(lbl)
        gi = find(strcmp(ch_labels, char(lbl)), 1);
    end

% ── Chunk dropdown management ─────────────────────────────────────────────
    function update_chunk_dropdown()
        npc  = n_per_chunk();
        nc   = n_chunks();
        strs = cell(nc,1);
        for ci = 1:nc
            s = (ci-1)*npc + 1;
            e = min(s+npc-1, n_chans);
            strs{ci} = sprintf('%s - %s', ch_labels{s}, ch_labels{e});
        end
        dd_chunk.String = strs;
        cur_chunk = max(1, min(cur_chunk, nc));
        dd_chunk.Value  = cur_chunk;
    end

    function update_chans_cb(~,~)
        v_rast = str2double(ed_raster.String);
        v_time = str2double(ed_timecourse.String);
        if ~isnan(v_rast) && v_rast > 0, cfg.n_chans_raster = round(v_rast); end
        if ~isnan(v_time) && v_time > 0, cfg.n_chans_timecourse = round(v_time); end
        custom_xlim = []; custom_ylim = [];
        update_chunk_dropdown();
        clip_sels_to_vis();
        update_plot();
    end

    function chunk_prev_cb(~,~)
        if dd_chunk.Value > 1
            dd_chunk.Value = dd_chunk.Value - 1;
            chunk_cb();
        end
    end

    function chunk_next_cb(~,~)
        if dd_chunk.Value < numel(dd_chunk.String)
            dd_chunk.Value = dd_chunk.Value + 1;
            chunk_cb();
        end
    end

% ── Trace Visibility & Scaling ────────────────────────────────────────────
    function toggle_vis_cb(src, ~)
        if strcmp(src.String, '✓')
            src.String = 'X';
            src.ForegroundColor = 'r';
        else
            src.String = '✓';
            src.ForegroundColor = 'k';
        end
        update_plot();
    end

    function scale_cb(type, delta)
        if strcmp(type, 'raw')
            raw_scale_val = max(0.1, raw_scale_val + delta);
            ed_raw_scale.String = num2str(raw_scale_val, '%.2f');
        elseif strcmp(type, 'clean')
            clean_scale_val = max(0.1, clean_scale_val + delta);
            ed_clean_scale.String = num2str(clean_scale_val, '%.2f');
        elseif strcmp(type, 'notch')
            notch_scale_val = max(0.1, notch_scale_val + delta);
            ed_notch_scale.String = num2str(notch_scale_val, '%.2f');
        end
        update_plot();
    end

% ── Axis Render Routines ──────────────────────────────────────────────────
    function update_plot()
        delete(findobj(fig,'Type','colorbar'));
        ax.Position = [LP_W+0.04, 0.09, 1-LP_W-0.06, 0.87];
        cla(ax);
        hold(ax,'on');
        vis   = get_vis_chans();
        n_vis = numel(vis);

        show_raw   = strcmp(btn_vis_raw.String, '✓');
        show_clean = strcmp(btn_vis_clean.String, '✓');
        show_notch = strcmp(btn_vis_notch.String, '✓');

        % ── RASTER MODE ──
        if strcmp(viewmode,'raster')
            set(lbl_cmap,'Visible','on');
            set(dd_cmap, 'Visible','on');
            set(txt_fs_display, 'Visible', 'off');
            
            vis_data = double(data_mat(vis,:));
            % Downsample Raster safely to prevent 1M+ column GPU crash
            n_samples = size(vis_data, 2);
            max_raster_pts = 4000;
            if n_samples > max_raster_pts
                ds_r = ceil(n_samples / max_raster_pts);
                vis_data = vis_data(:, 1:ds_r:end);
            end
            
            for i = 1:size(vis_data,1)
                s = std(vis_data(i,:));
                if s < eps('single'), s = 1; end
                vis_data(i,:) = (vis_data(i,:) - mean(vis_data(i,:))) / s;
            end
            
            % Map raster spatially to X limit range instead of 1M+ element discrete array
            x_lims = [time_vec(1), time_vec(end)];
            imagesc(ax, x_lims, [1 n_vis], vis_data);
            colormap(ax, cur_cmap);
            
            try 
                cbh = colorbar(ax); 
                cbh.FontSize = 8; 
                cbh.Ruler.TickLabelFormat = '%.1f'; 
                cbh.Ruler.Exponent = 0; 
            catch
            end
            
            for ai = 1:height(artifact)
                gi = lbl2idx(artifact.label(ai));
                if isempty(gi), continue; end
                li = find(vis==gi,1);
                if isempty(li), continue; end
                rect_outline(artifact.starts(ai), artifact.ends(ai), ...
                             li-0.5, li+0.5, [0.85 0.07 0.07], 2.5);
            end
            for si = 1:numel(SEL)
                [li1,li2] = sel2local(SEL(si),vis,n_vis);
                if isnan(li1), continue; end
                rect_outline(SEL(si).t_start, SEL(si).t_end, ...
                             li1-0.5, li2+0.5, [0.25 0.25 0.25], 2.0);
            end

        % ── TIMECOURSE MODE ──
        else
            set(lbl_cmap,'Visible','off');
            set(dd_cmap, 'Visible','off');
            set(txt_fs_display, 'Visible', 'on');
            
            % Compute original sample rate safely
            n_pts = numel(time_vec);
            dt = mean(diff(time_vec(1:min(1000, end))));
            orig_fs = 1 / dt;
            
            % Downsample if rate exceeds the configured maximum limit
            if orig_fs > op.max_timecourse_sample_rate
                ds = max(1, round(orig_fs / op.max_timecourse_sample_rate));
                idx_ = 1:ds:n_pts;
                display_fs = orig_fs / ds;
            else
                idx_ = 1:n_pts;
                display_fs = orig_fs;
            end
            
            % Force mathematically strictly to 1xM row vectors
            t_d = time_vec(idx_);
            t_d = t_d(:)'; 
            
            set(txt_fs_display, 'String', sprintf('Sample Rate: %.1f Hz', display_fs));
            
            for ai = 1:height(artifact)
                gi = lbl2idx(artifact.label(ai));
                if isempty(gi), continue; end
                li = find(vis==gi,1);
                if isempty(li), continue; end
                rect_fill(artifact.starts(ai), artifact.ends(ai), ...
                          li-0.48, li+0.48, [1 0.67 0.67], [0.80 0.05 0.05], 1.8);
            end
            for si = 1:numel(SEL)
                [li1,li2] = sel2local(SEL(si),vis,n_vis);
                if isnan(li1), continue; end
                rect_fill(SEL(si).t_start, SEL(si).t_end, ...
                          li1-0.48, li2+0.48, [0.78 0.78 0.78], [0.33 0.33 0.33], 1.5);
            end
            for ki = 1:n_vis-1
                yline(ax, ki+0.5, 'Color',[0.79 0.79 0.79], 'LineWidth',0.5);
            end
            
            for ki = 1:n_vis
                
                % Process and plot Raw trace dynamically
                if show_raw && ~isempty(data_mat)
                    len_raw = size(data_mat, 2);
                    idx_raw = idx_(idx_ <= len_raw);
                    if ~isempty(idx_raw)
                        sig = double(data_mat(vis(ki), idx_raw));
                        sig = sig(:)'; % Force row vector
                        s = std(sig);
                        if s < eps('single'), s = 1; end
                        sig_n = (sig - mean(sig)) ./ (6*s) * raw_scale_val + ki;
                        plot(ax, t_d(1:numel(idx_raw)), sig_n, 'Color','k', 'LineWidth',0.5);
                    end
                end
                
                % Process and plot Clean trace dynamically
                if show_clean && ~isempty(data_mat_clean)
                    len_clean = size(data_mat_clean, 2);
                    idx_clean = idx_(idx_ <= len_clean);
                    if ~isempty(idx_clean)
                        sig_clean = double(data_mat_clean(vis(ki), idx_clean));
                        sig_clean = sig_clean(:)'; % Force row vector
                        s_clean = std(sig_clean);
                        if s_clean < eps('single'), s_clean = 1; end
                        sig_clean_n = (sig_clean - mean(sig_clean)) ./ (6*s_clean) * clean_scale_val + ki;
                        plot(ax, t_d(1:numel(idx_clean)), sig_clean_n, 'Color','b', 'LineWidth',0.5);
                    end
                end
                
                % Process and plot Notch trace dynamically
                if show_notch && ~isempty(data_mat_notch)
                    len_notch = size(data_mat_notch, 2);
                    idx_notch = idx_(idx_ <= len_notch);
                    if ~isempty(idx_notch)
                        sig_notch = double(data_mat_notch(vis(ki), idx_notch));
                        sig_notch = sig_notch(:)'; % Force row vector
                        s_notch = std(sig_notch);
                        if s_notch < eps('single'), s_notch = 1; end
                        sig_notch_n = (sig_notch - mean(sig_notch)) ./ (6*s_notch) * notch_scale_val + ki;
                        plot(ax, t_d(1:numel(idx_notch)), sig_notch_n, 'Color',[0.6 0 0], 'LineWidth',0.5);
                    end
                end
                
            end
        end

        % ── Shared Aesthetics Configuration ──
        yticks(ax, 1:n_vis);
        yticklabels(ax, ch_labels(vis));
        ax.TickLabelInterpreter = 'none';
        
        if isempty(custom_ylim)
            ylim(ax, [0.5, n_vis+0.5]);
        else
            ylim(ax, custom_ylim);
        end
        
        if isempty(custom_xlim)
            xlim(ax, [time_vec(1), time_vec(end)]);
        else
            xlim(ax, custom_xlim);
        end
        
        xlabel(ax,'Time (s)');
        ylabel(ax,'Channel');
        ax.YDir    = 'reverse';   
        ax.FontSize = 8;
        ax.Layer   = 'top';
        
        ax.XAxis.TickLabelFormat = '%.1f';
        ax.XAxis.Exponent = 0;
        
        drawnow limitrate;
        update_selection_info();
    end

    function update_selection_info()
        if isempty(SEL)
            txt_sel_info.String = 'Selection: None';
        else
            dur = SEL(1).t_end - SEL(1).t_start;
            n_ch = abs(SEL(1).ch_end - SEL(1).ch_start) + 1;
            txt_sel_info.String = sprintf('Selection: Duration = %.3f s, Channels selected = %d', dur, n_ch);
        end
    end

    function rect_outline(t0,t1,y0,y1,clr,lw)
        plot(ax,[t0 t1 t1 t0 t0],[y0 y0 y1 y1 y0], ...
            '-','Color',clr,'LineWidth',lw,'HitTest','off');
    end

    function rect_fill(t0,t1,y0,y1,fc,ec,lw)
        patch(ax,[t0 t1 t1 t0],[y0 y0 y1 y1], ...
            fc,'EdgeColor',ec,'LineWidth',lw,'FaceAlpha',0.50,'HitTest','off');
    end

    function [li1,li2] = sel2local(sel,vis,n_vis)
        inr = vis(vis >= sel.ch_start & vis <= sel.ch_end);
        if isempty(inr), li1=NaN; li2=NaN; return; end
        li1 = find(vis==inr(1),  1);
        li2 = find(vis==inr(end),1);
        li1 = max(1,min(n_vis,li1));
        li2 = max(1,min(n_vis,li2));
    end

% ── Hotkey Commands / Window Callbacks ────────────────────────────────────
    function key_press_cb(~, event)
        switch event.Character
            case 'a'
                add_artifact_cb();
            case 'r'
                remove_artifact_cb();
            case 'q'
                zoom_sel_cb();
            case 'w'
                zoom_full_cb();
            case 'm'
                viewmode_cb();
            case 'o'
                load_artifact_cb();
            case 'k'
                save_artifact_cb();
        end
    end

    function zoom_sel_cb(~,~)
        if isempty(SEL), return; end
        vis = get_vis_chans();
        n_vis = numel(vis);
        [li1, li2] = sel2local(SEL(1), vis, n_vis);
        if isnan(li1), return; end
        custom_xlim = [SEL(1).t_start, SEL(1).t_end];
        custom_ylim = [li1 - 0.5, li2 + 0.5];
        update_plot();
    end

    function zoom_full_cb(~,~)
        custom_xlim = [];
        custom_ylim = [];
        update_plot();
    end

    function update_viewmode_button()
        % Uses HTML to independently format both strings on the same button
        if strcmp(viewmode, 'timecourse')
            str = '<html><font color="#008000" size="5"><b>TIMECOURSE</b></font><font color="#000000" size="3"> // raster (m)</font></html>';
        else
            str = '<html><font color="#000000" size="3">timecourse // </font><font color="#008000" size="5"><b>RASTER (m)</b></font></html>';
        end
        set(btn_view, 'String', str);
    end

    function viewmode_cb(~, ~)
        if strcmp(viewmode, 'timecourse')
            target_mode = 'raster';
        else
            target_mode = 'timecourse';
        end
        
        old_vis = get_vis_chans();
        viewmode = target_mode;
        
        update_viewmode_button();
        custom_xlim = []; custom_ylim = [];
        cur_chunk = best_chunk_for(old_vis);
        update_chunk_dropdown();
        clip_sels_to_vis();
        update_plot();
    end

    function best = best_chunk_for(old_vis)
        npc=n_per_chunk(); nc=n_chunks();
        best=1; bn=-1;
        for ci=1:nc
            s=(ci-1)*npc+1; e=min(s+npc-1,n_chans);
            ov=numel(intersect(old_vis,s:e));
            if ov>bn, bn=ov; best=ci; end
        end
    end

    function clip_sels_to_vis()
        if isempty(SEL), return; end
        vis=get_vis_chans();
        keep=true(1,numel(SEL));
        for si=1:numel(SEL)
            inr=vis(vis>=SEL(si).ch_start & vis<=SEL(si).ch_end);
            if isempty(inr), keep(si)=false;
            else, SEL(si).ch_start=inr(1); SEL(si).ch_end=inr(end); end
        end
        SEL=SEL(keep);
    end

    function chunk_cb(~,~)
        cur_chunk=dd_chunk.Value;
        custom_xlim = []; custom_ylim = [];
        update_plot();
    end

    function cmap_cb(~,~)
        cur_cmap=CMAPS{dd_cmap.Value};
        update_plot();
    end

% ── Mouse callbacks ───────────────────────────────────────────────────────
    function btn_down_cb(~,~)
        cp=ax_coords_clamped();
        if isempty(cp), return; end
        is_dragging=true;
        drag_x0=cp(1); drag_y0=cp(2);
        if isgraphics(rband_h), delete(rband_h); end
        rband_h = patch(ax, ...
            'XData', repmat(drag_x0,1,4), ...
            'YData', repmat(drag_y0,1,4), ...
            'FaceColor','none','EdgeColor',[0.20 0.20 0.20], ...
            'LineWidth',1.5,'LineStyle','--','HitTest','off');
    end

    function btn_motion_cb(~,~)
        if ~is_dragging, return; end
        cp=ax_coords_raw();
        if isempty(cp), return; end
        x0=min(drag_x0,cp(1)); x1=max(drag_x0,cp(1));
        y0=min(drag_y0,cp(2)); y1=max(drag_y0,cp(2));
        if isgraphics(rband_h)
            set(rband_h,'XData',[x0 x1 x1 x0],'YData',[y0 y0 y1 y1]);
        end
    end

    function btn_up_cb(~,~)
        if ~is_dragging, return; end
        is_dragging=false;
        if isgraphics(rband_h), delete(rband_h); end
        rband_h=[];
        cp=ax_coords_raw();
        if isempty(cp), return; end
        t_lo=max(min(drag_x0,cp(1)), time_vec(1));
        t_hi=min(max(drag_x0,cp(1)), time_vec(end));
        if t_hi<=t_lo, return; end
        y_lo=min(drag_y0,cp(2));
        y_hi=max(drag_y0,cp(2));
        vis=get_vis_chans(); n_vis=numel(vis);
        li_lo=max(1,     ceil(y_lo - 0.5));
        li_hi=min(n_vis, floor(y_hi + 0.5));
        if li_lo>li_hi
            li_lo=max(1,min(n_vis,round((y_lo+y_hi)/2)));
            li_hi=li_lo;
        end
        SEL = struct('t_start',t_lo,'t_end',t_hi, ...
                     'ch_start',vis(li_lo),'ch_end',vis(li_hi));
        update_plot();
    end

    function cp=ax_coords_clamped()
        cp=ax_coords_raw();
        if isempty(cp), return; end
        xl=sort(xlim(ax)); yl=sort(ylim(ax));
        if cp(1)<xl(1)||cp(1)>xl(2)||cp(2)<yl(1)||cp(2)>yl(2), cp=[]; end
    end

    function cp=ax_coords_raw()
        pt=ax.CurrentPoint;
        if isempty(pt), cp=[]; return; end
        cp=[pt(1,1), pt(1,2)];
    end

% ── Add artifact ──────────────────────────────────────────────────────────
    function add_artifact_cb(~,~)
        if isempty(SEL)
            msgbox('No region selected.','Add Artifact','help'); return;
        end
        rows=sels_to_rows(SEL);
        SEL=struct('t_start',{},'t_end',{},'ch_start',{},'ch_end',{});
        if height(rows)==0, update_plot(); return; end
        combine_artifact_tables(rows);
    end

    function rows=sels_to_rows(sels)
        rows=mk_empty_table();
        for si=1:numel(sels)
            s=sels(si);
            for gi=s.ch_start:s.ch_end
                rows=[rows; mk_row(0,s.t_start,s.t_end,ch_labels{gi})]; 
            end
        end
    end

% ── Remove artifact ───────────────────────────────────────────────────────
    function remove_artifact_cb(~,~)
        if isempty(SEL)
            msgbox('No region selected.','Remove Artifact','help'); return;
        end
        if height(artifact)==0
            msgbox('Artifact table is empty.','Remove Artifact','help'); return;
        end
        aff=false(height(artifact),1);
        for ai=1:height(artifact)
            gi=lbl2idx(artifact.label(ai));
            if isempty(gi), continue; end
            for si=1:numel(SEL)
                sel=SEL(si);
                if gi<sel.ch_start||gi>sel.ch_end, continue; end
                if artifact.ends(ai)<=sel.t_start, continue; end
                if artifact.starts(ai)>=sel.t_end, continue; end
                aff(ai)=true; break;
            end
        end
        if ~any(aff)
            msgbox('No overlap with existing artifacts.','Remove Artifact','help'); return;
        end
        aff_idx=find(aff);
        if numel(aff_idx)>1 && is_multi_group(aff_idx)
            btn=questdlg('Remove multiple groups of time x channels?', ...
                'Confirm Removal','Yes','No','No');
            if isempty(btn)||strcmp(btn,'No'), return; end
        end
        new_art=mk_empty_table();
        for ai=1:height(artifact)
            gi=lbl2idx(artifact.label(ai));
            if ~aff(ai)||isempty(gi)
                new_art=[new_art; artifact(ai,:)]; 
                continue;
            end
            ivs=[artifact.starts(ai), artifact.ends(ai)];
            for si=1:numel(SEL)
                sel=SEL(si);
                if gi<sel.ch_start||gi>sel.ch_end, continue; end
                ivs=sub_interval(ivs,sel.t_start,sel.t_end);
                if isempty(ivs), break; end
            end
            for ri=1:size(ivs,1)
                if ivs(ri,2)-ivs(ri,1)>1e-10
                    new_art=[new_art; mk_row(0,ivs(ri,1),ivs(ri,2),ch_labels{gi})]; 
                end
            end
        end
        artifact=sortrows(new_art,{'starts','label'});
        artifact.id=(1:height(artifact))';
        update_plot();
    end

    function out=sub_interval(ivs,rm_s,rm_e)
        out=zeros(0,2);
        for k=1:size(ivs,1)
            a=ivs(k,1); b=ivs(k,2);
            if     rm_e<=a||rm_s>=b,       out(end+1,:)=[a,b];        
            elseif rm_s<=a&&rm_e>=b,       
            elseif rm_s<=a,                out(end+1,:)=[rm_e,b];     
            elseif rm_e>=b,                out(end+1,:)=[a,rm_s];     
            else, out(end+1,:)=[a,rm_s]; out(end+1,:)=[rm_e,b];       
            end
        end
    end

    function result=is_multi_group(idx)
        n=numel(idx);
        if n<=1, result=false; return; end
        adj=false(n);
        for ii=1:n
            gi_i=lbl2idx(artifact.label(idx(ii)));
            for jj=ii+1:n
                gi_j=lbl2idx(artifact.label(idx(jj)));
                ch_ok=~isempty(gi_i)&&~isempty(gi_j)&&abs(gi_i-gi_j)<=1;
                t_ok = artifact.starts(idx(ii)) < artifact.ends(idx(jj)) && ...
                       artifact.starts(idx(jj)) < artifact.ends(idx(ii));
                if ch_ok&&t_ok, adj(ii,jj)=true; adj(jj,ii)=true; end
            end
        end
        visited=false(n,1); n_comp=0;
        for k=1:n
            if ~visited(k)
                n_comp=n_comp+1;
                if n_comp>1, result=true; return; end
                q=k;
                while ~isempty(q)
                    cur=q(1); q(1)=[];
                    if visited(cur), continue; end
                    visited(cur)=true;
                    q=[q, find(adj(cur,:)&~visited')]; 
                end
            end
        end
        result=n_comp>1;
    end

    function combine_artifact_tables(new_rows)
        artifact=[artifact; new_rows];
        u_lbls=unique(artifact.label);
        merged=mk_empty_table();
        for ui=1:numel(u_lbls)
            lbl  =u_lbls(ui);
            sub_t=sortrows(artifact(artifact.label==lbl,:),'starts');
            ms=sub_t.starts(1); me=sub_t.ends(1);
            for ri=2:height(sub_t)
                if sub_t.starts(ri)<=me+1e-10
                    me=max(me,sub_t.ends(ri));
                else
                    merged=[merged; mk_row(0,ms,me,char(lbl))]; 
                    ms=sub_t.starts(ri); me=sub_t.ends(ri);
                end
            end
            merged=[merged; mk_row(0,ms,me,char(lbl))]; 
        end
        merged   =sortrows(merged,{'starts','label'});
        merged.id=(1:height(merged))';
        artifact =merged;
        update_plot();
    end

% ── Load Artifact File Master ─────────────────────────────────────────────
    function load_artifact_file(fpath)
        try
            al = readtable(fpath, 'FileType','text','Delimiter','\t','TextType','string');
            probs = validate_art_table(al);
            if ~isempty(probs)
                errordlg(['Validation failed:' newline strjoin(probs,newline)], ...
                    'Invalid Artifact Table');
                return;
            end
            req={'starts','ends','duration','label'};
            xtra=setdiff(al.Properties.VariableNames,req);
            for k=1:numel(xtra), al.(xtra{k})=[]; end
            al.label=string(al.label);
            al.id   =zeros(height(al),1,'double');
            al=al(:,{'id','starts','ends','duration','label'});
            combine_artifact_tables(al);
        catch ME
            errordlg(sprintf('Cannot read file:\n%s',ME.message),'Load Error');
        end
    end

% ── Manual Load Artifact Mask ─────────────────────────────────────────────
    function load_artifact_cb(~,~)
        if height(artifact)>0
            btn=questdlg( ...
                'Artifact table is not empty. Add loaded table to current artifacts?', ...
                'Load Artifact Mask','Yes','No','Yes');
            if isempty(btn)||strcmp(btn,'No'), return; end
        end
        annot_dir=fullfile('Y:\DBS','derivatives',['sub-' cfg.sub],'annot');
        [fname,fdir]=uigetfile( ...
            {'*artifact*.tsv','Artifact TSV (*artifact*.tsv)'; ...
             '*.tsv',         'All TSV (*.tsv)'}, ...
            'Load Artifact Mask', annot_dir);
        if isequal(fname,0), return; end
        
        tok = regexp(fname, '^sub-([^_]+)_', 'tokens');
        if ~isempty(tok)
            file_sub = tok{1}{1};
            if ~strcmp(file_sub, cfg.sub)
                ans_btn = questdlg(sprintf('Subject mismatch!\nLoaded file subject: %s\nCurrent workspace subject: %s\n\nProceed anyway?', file_sub, cfg.sub), 'Subject Mismatch', 'Proceed', 'Cancel', 'Cancel');
                if isempty(ans_btn) || strcmp(ans_btn, 'Cancel')
                    return;
                end
            end
        end

        load_artifact_file(fullfile(fdir, fname));
    end

    function probs=validate_art_table(T)
        probs={};
        if height(T)==0
            probs{end+1}='  * Table has zero rows.'; return;
        end
        req_cols={'starts','ends','duration','label'};
        for ci=1:numel(req_cols)
            col=req_cols{ci};
            if ~ismember(col,T.Properties.VariableNames)
                probs{end+1}=sprintf('  * Missing column: ''%s''.',col);
            else
                v=T.(col);
                if ismember(col,{'starts','ends','duration'})&&~isnumeric(v)
                    probs{end+1}=sprintf('  * ''%s'' must be numeric (got %s).',col,class(v));
                end
                if strcmp(col,'label')&&~(isstring(v)||iscell(v))
                    probs{end+1}=sprintf('  * ''label'' must be string/cell (got %s).',class(v));
                end
            end
        end
        if ~isempty(probs), return; end
        if any(abs((T.ends-T.starts)-T.duration)>1e-9)
            probs{end+1}='  * duration != ends-starts for some rows.';
        end
        if any(T.duration<=0|isnan(T.duration))
            probs{end+1}='  * Some durations are <= 0 or NaN.';
        end
        ls=string(T.label);
        bad_lbl=~ismember(ls,string(ch_labels));
        if any(bad_lbl)
            ul=unique(ls(bad_lbl));
            probs{end+1}=sprintf('  * %d label(s) not in FieldTrip object: %s', ...
                numel(ul),strjoin(ul,', '));
        end
    end

% ── Save Artifact Mask ────────────────────────────────────────────────────
    function save_artifact_cb(~,~)
        annot_dir   =fullfile('Y:\DBS','derivatives',['sub-' cfg.sub],'annot');
        default_name=sprintf('sub-%s_ses-intraop_task-%s_artifact-manual.tsv',cfg.sub,cfg.task);
        [fname,fdir]=uiputfile({'*.tsv','Tab-separated values (*.tsv)'}, ...
            'Save Artifact Mask',fullfile(annot_dir,default_name));
        if isequal(fname,0), return; end
        try
            writetable(artifact,fullfile(fdir,fname), ...
                'FileType','text','Delimiter','\t');
            msgbox(sprintf('Saved:\n%s',fullfile(fdir,fname)),'Saved');
        catch ME
            errordlg(sprintf('Save failed:\n%s',ME.message),'Save Error');
        end
    end

% ── Data Table Formatting Helpers ─────────────────────────────────────────
    function T=mk_empty_table()
        T=table(zeros(0,1,'double'),zeros(0,1,'double'), ...
                zeros(0,1,'double'),zeros(0,1,'double'),strings(0,1), ...
                'VariableNames',{'id','starts','ends','duration','label'});
    end

    function row=mk_row(id_,t0_,t1_,lbl_char)
        row=table(double(id_),double(t0_),double(t1_),double(t1_-t0_), ...
                  string(lbl_char), ...
                  'VariableNames',{'id','starts','ends','duration','label'});
    end

end