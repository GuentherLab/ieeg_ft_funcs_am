function annotate_manual_ephys_artifacts(op)
%ANNOTATE_MANUAL_EPHYS_ARTIFACTS  Interactive GUI for manual annotation of
%   intra-operative ephys artifacts in a FieldTrip raw structure.
%
% Usage
%   annotate_manual_ephys_artifacts()
%   annotate_manual_ephys_artifacts(op)
%
% op fields (all optional)
%   .sub                – Subject ID string, e.g. 'DM1005'  (prompted if absent)
%   .n_chans_raster     – Channels shown per screen in raster view     [40]
%   .n_chans_timecourse – Channels shown per screen in timecourse view [20]
%
% Input file
%   Y:\DBS\derivatives\sub-<SUB>\fieldtrip\sub-<SUB>_ses-intraop_task-smsl_ft-raw.mat
%   Must contain a FieldTrip raw struct with fields: trial, time, label.
%   data_mat = D.trial{1}  [nCh x nTime]
%   time_vec = D.time{1}   [1   x nTime]  (global time coordinates, seconds)
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
if nargin < 1 || isempty(op), op = struct(); end
if ~isfield(op,'n_chans_raster'),      op.n_chans_raster     = 40; end
if ~isfield(op,'n_chans_timecourse'),  op.n_chans_timecourse = 20; end
%==========================================================================
%% 1.  Subject ID
%==========================================================================
if ~isfield(op,'sub') || isempty(op.sub)
    ans_ = inputdlg('Enter subject ID:','Subject ID',1,{'DM10'});
    if isempty(ans_), return; end
    op.sub = strtrim(ans_{1});
end
%==========================================================================
%% 2.  Load FieldTrip file
%==========================================================================
ft_path = fullfile('Y:\DBS','derivatives', ...
    ['sub-' op.sub], 'fieldtrip', ...
    ['sub-' op.sub '_ses-intraop_task-smsl_ft-raw.mat']);
fprintf('\n=== FieldTrip file ===\n');
fprintf('  Path : %s\n', ft_path);
dinfo = dir(ft_path);
if isempty(dinfo)
    errordlg(sprintf('File not found:\n%s',ft_path),'File Not Found');
    return;
end
fprintf('  Size : %.2f MB\n', dinfo.bytes/1e6);
t_load0 = datetime('now');
fprintf('  Load started : %s\n', char(t_load0,'HH:mm:ss.SSS'));
loaded = load(ft_path);
t_load1 = datetime('now');
fprintf('  Load done    : %s  (%.2f s)\n\n', ...
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
    errordlg('No FieldTrip raw struct (trial/label/time) found in file.','Load Error');
    return;
end
data_mat  = D.trial{1};         % [nCh x nTime]
time_vec  = D.time{1};          % [1   x nTime]
ch_labels = cellstr(D.label);   % {nCh x 1}
n_chans   = size(data_mat,1);
%==========================================================================
%% 3.  Shared state  (accessed/modified by nested functions)
%==========================================================================
viewmode  = 'timecourse';    % 'raster' | 'timecourse' (Default changed)
cur_chunk = 1;
% Custom zoom limits
custom_xlim = [];
custom_ylim = [];
% Artifact table
artifact = table( ...
    zeros(0,1,'double'), zeros(0,1,'double'), ...
    zeros(0,1,'double'), zeros(0,1,'double'), strings(0,1), ...
    'VariableNames',{'id','starts','ends','duration','label'});
% Current gray selections – struct array
%   .t_start, .t_end   : global time [s]
%   .ch_start, .ch_end : global channel index (1-based)
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
MAX_TC_PTS = 5000;  % max displayed time-points per trace (timecourse mode)
%==========================================================================
%% 4.  Build figure & controls
%==========================================================================
LP_W = 0.17;   % fractional figure width for left panel
fig = figure( ...
    'Name',          ['Artifact Annotator  |  sub-' op.sub], ...
    'NumberTitle',   'off', ...
    'Units',         'normalized', ...
    'OuterPosition', [0.01 0.03 0.98 0.94], ...
    'Color',         [0.93 0.93 0.93], ...
    'WindowButtonDownFcn',   @btn_down_cb, ...
    'WindowButtonMotionFcn', @btn_motion_cb, ...
    'WindowButtonUpFcn',     @btn_up_cb, ...
    'WindowKeyPressFcn',     @key_press_cb);
% ── Left panel ────────────────────────────────────────────────────────────
lp = uipanel(fig, ...
    'Title','Controls','Units','normalized', ...
    'Position',[0 0 LP_W 1], ...
    'BackgroundColor',[0.87 0.87 0.87],'FontSize',9);
% Stacking helpers
ny_ = 0.985;  ch_ = 0.035;  dh_ = 0.040;

% View mode Toggle
yp_ = ny_-ch_; ny_ = ny_-dh_;
btn_viewmode = uicontrol(lp,'Style','pushbutton', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@viewmode_cb);

% Colormap
yp_ = ny_-ch_; ny_ = ny_-dh_;
lbl_cmap = uicontrol(lp,'Style','text','String','Colormap:', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
yp_ = ny_-ch_; ny_ = ny_-dh_;
dd_cmap = uicontrol(lp,'Style','popupmenu','String',CMAPS, ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@cmap_cb);

% Channel count modifiers
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','N Visible Chans Raster:', ...
    'Units','normalized','Position',[0.05 yp_ 0.50 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
ed_raster = uicontrol(lp,'Style','edit','String',num2str(op.n_chans_raster), ...
    'Units','normalized','Position',[0.60 yp_ 0.35 ch_]);
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','text','String','N Visible Chans Timecourse:', ...
    'Units','normalized','Position',[0.05 yp_ 0.50 ch_], ...
    'HorizontalAlignment','left','BackgroundColor',[0.87 0.87 0.87]);
ed_timecourse = uicontrol(lp,'Style','edit','String',num2str(op.n_chans_timecourse), ...
    'Units','normalized','Position',[0.60 yp_ 0.35 ch_]);
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Update Chans', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@update_chans_cb);
ny_ = ny_ - 0.010;   % spacer

% Channel group
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

% Artifact action buttons
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Add artifact (a)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_], ...
    'FontWeight','bold','Callback',@add_artifact_cb);
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Remove selected artifact (r)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@remove_artifact_cb);
ny_ = ny_ - 0.018;   % spacer

% File buttons
yp_ = ny_-ch_; ny_ = ny_-dh_;
uicontrol(lp,'Style','pushbutton','String','Load artifact mask (o)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@load_artifact_cb);
yp_ = ny_-ch_; % last row
uicontrol(lp,'Style','pushbutton','String','Save artifact mask (k)', ...
    'Units','normalized','Position',[0.05 yp_ 0.90 ch_],'Callback',@save_artifact_cb);

% ── Main axes ─────────────────────────────────────────────────────────────
ax = axes(fig,'Units','normalized', ...
    'Position',[LP_W+0.04, 0.08, 1-LP_W-0.06, 0.88]);
%==========================================================================
%% 5.  Initialise display
%==========================================================================
update_viewmode_button();
update_chunk_dropdown();
update_plot();

%==========================================================================
%%   ═══════════  N E S T E D   F U N C T I O N S  ═══════════
%==========================================================================

% ── Simple helpers ────────────────────────────────────────────────────────
    function n = n_per_chunk()
        if strcmp(viewmode,'raster'), n = op.n_chans_raster;
        else,                         n = op.n_chans_timecourse; end
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

% ── Chunk dropdown ────────────────────────────────────────────────────────
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
        if ~isnan(v_rast) && v_rast > 0, op.n_chans_raster = round(v_rast); end
        if ~isnan(v_time) && v_time > 0, op.n_chans_timecourse = round(v_time); end
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

% ── Main plot ─────────────────────────────────────────────────────────────
    function update_plot()
        delete(findobj(fig,'Type','colorbar'));
        ax.Position = [LP_W+0.04, 0.08, 1-LP_W-0.06, 0.88];
        cla(ax);
        hold(ax,'on');
        vis   = get_vis_chans();
        n_vis = numel(vis);

        % ── RASTER ───────────────────────────────────────────────────────
        if strcmp(viewmode,'raster')
            set(lbl_cmap,'Visible','on');
            set(dd_cmap, 'Visible','on');
            % Standardise each row relative to itself to avoid blowouts
            vis_data = double(data_mat(vis,:));
            for i = 1:size(vis_data,1)
                s = std(vis_data(i,:));
                if s < eps('single'), s = 1; end
                vis_data(i,:) = (vis_data(i,:) - mean(vis_data(i,:))) / s;
            end
            imagesc(ax, time_vec, 1:n_vis, vis_data);
            colormap(ax, cur_cmap);
            try cbh = colorbar(ax); cbh.FontSize = 8; catch; end
            
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

        % ── TIMECOURSE ───────────────────────────────────────────────────
        else
            set(lbl_cmap,'Visible','off');
            set(dd_cmap, 'Visible','off');
            
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
            nt = numel(time_vec);
            if nt > MAX_TC_PTS
                ds   = floor(nt / MAX_TC_PTS);
                t_d  = time_vec(1:ds:end);
                idx_ = 1:ds:nt;
            else
                t_d  = time_vec;
                idx_ = 1:nt;
            end
            for ki = 1:n_vis
                sig   = double(data_mat(vis(ki), idx_));
                s     = std(sig(:));
                if s < eps('single'), s = 1; end
                sig_n = (sig - mean(sig)) ./ (6*s) + ki;
                plot(ax, t_d, sig_n, 'Color','k', 'LineWidth',0.5);
            end
        end

        % ── Common formatting ─────────────────────────────────────────────
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
        ax.YDir    = 'reverse';   % channel 1 at top
        ax.FontSize = 8;
        ax.Layer   = 'top';
        drawnow limitrate;
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

% ── Key Press / Zoom Callbacks ────────────────────────────────────────────
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

% ── View-mode / chunk / colormap callbacks ────────────────────────────────
    function update_viewmode_button()
        % Uses HTML to give the active state a background tint, green text, and larger font size.
        if strcmp(viewmode, 'timecourse')
            btn_str = '<html><span style="background-color:#d4edda; color:#007A33; font-size:1.15em;"><b>&nbsp;Timecourse&nbsp;</b></span> // Raster (m)</html>';
        else
            btn_str = '<html>Timecourse // <span style="background-color:#d4edda; color:#007A33; font-size:1.15em;"><b>&nbsp;Raster&nbsp;</b></span> (m)</html>';
        end
        set(btn_viewmode, 'String', btn_str);
    end

    function viewmode_cb(~,~)
        old_vis = get_vis_chans();
        if strcmp(viewmode, 'timecourse')
            viewmode = 'raster';
        else
            viewmode = 'timecourse';
        end
        
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
        % Replaces prior selection entirely rather than appending
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
                rows=[rows; mk_row(0,s.t_start,s.t_end,ch_labels{gi})]; %#ok<AGROW>
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
                new_art=[new_art; artifact(ai,:)]; %#ok<AGROW>
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
                    new_art=[new_art; mk_row(0,ivs(ri,1),ivs(ri,2),ch_labels{gi})]; %#ok<AGROW>
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
            if     rm_e<=a||rm_s>=b,       out(end+1,:)=[a,b];        %#ok<AGROW>
            elseif rm_s<=a&&rm_e>=b,       
            elseif rm_s<=a,                out(end+1,:)=[rm_e,b];     %#ok<AGROW>
            elseif rm_e>=b,                out(end+1,:)=[a,rm_s];     %#ok<AGROW>
            else, out(end+1,:)=[a,rm_s]; out(end+1,:)=[rm_e,b];       %#ok<AGROW>
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
                    q=[q, find(adj(cur,:)&~visited')]; %#ok<AGROW>
                end
            end
        end
        result=n_comp>1;
    end

% ── combine_artifact_tables ───────────────────────────────────────────────
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
                    merged=[merged; mk_row(0,ms,me,char(lbl))]; %#ok<AGROW>
                    ms=sub_t.starts(ri); me=sub_t.ends(ri);
                end
            end
            merged=[merged; mk_row(0,ms,me,char(lbl))]; %#ok<AGROW>
        end
        merged   =sortrows(merged,{'starts','label'});
        merged.id=(1:height(merged))';
        artifact =merged;
        update_plot();
    end

% ── Load artifact mask ────────────────────────────────────────────────────
    function load_artifact_cb(~,~)
        if height(artifact)>0
            btn=questdlg( ...
                'Artifact table is not empty. Add loaded table to current artifacts?', ...
                'Load Artifact Mask','Yes','No','Yes');
            if isempty(btn)||strcmp(btn,'No'), return; end
        end
        annot_dir=fullfile('Y:\DBS','derivatives',['sub-' op.sub],'annot');
        [fname,fdir]=uigetfile( ...
            {'*artifact*.tsv','Artifact TSV (*artifact*.tsv)'; ...
             '*.tsv',         'All TSV (*.tsv)'}, ...
            'Load Artifact Mask', annot_dir);
        if isequal(fname,0), return; end
        
        % Subject mismatch check
        tok = regexp(fname, '^sub-([^_]+)_', 'tokens');
        if ~isempty(tok)
            file_sub = tok{1}{1};
            if ~strcmp(file_sub, op.sub)
                ans_btn = questdlg(sprintf('Subject mismatch!\nLoaded file subject: %s\nCurrent workspace subject: %s\n\nProceed anyway?', file_sub, op.sub), 'Subject Mismatch', 'Proceed', 'Cancel', 'Cancel');
                if isempty(ans_btn) || strcmp(ans_btn, 'Cancel')
                    return;
                end
            end
        end

        try
            al=readtable(fullfile(fdir,fname), ...
                'FileType','text','Delimiter','\t','TextType','string');
        catch ME
            errordlg(sprintf('Cannot read file:\n%s',ME.message),'Load Error');
            return;
        end
        probs=validate_art_table(al);
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
%         bad_t=T.starts<time_vec(1)|T.ends>time_vec(end);
%         if any(bad_t)
%             probs{end+1}=sprintf('  * %d row(s) outside time range [%.4f, %.4f] s.', ...
%                 sum(bad_t),time_vec(1),time_vec(end));
%         end
    end

% ── Save artifact mask ────────────────────────────────────────────────────
    function save_artifact_cb(~,~)
        annot_dir   =fullfile('Y:\DBS','derivatives',['sub-' op.sub],'annot');
        default_name=sprintf('sub-%s_ses-intraop_task-smsl_artifact-manual.tsv',op.sub);
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

% ── Table utilities ───────────────────────────────────────────────────────
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

%% claude/gemini prompts:
% create a matlab function ‘annotate_manual_ephys_artifacts’ that does the following.
% i have a fieldtrip variable stored in a file with filepath formatted as: 
% Y:\DBS\derivatives\sub-[SUBJECT]\fieldtrip\sub-[SUBJECT]_ses-intraop_task-smsl_ft-raw_trial.mat
% 
% function should take ‘op’ structure. if op.sub ([SUBJECT] above)  is not supplied, make a popup, prefilled with ‘DM1005’, where user can fill it in. 
% 
% open the fieldtrip file.within D.trial{1}, rows are electrodes, columns are timepoints. look up file formatting in fieldtrip to understand. before loading, display the filepath of the fieldtirp file and its size, and the time loading started. Once it’s loaded, display the time it took to load it. 
% 
% when it’s loaded, create a gui which displays the timecourses of a set of channels. there should be 2 main viewmodes: timecourse and raster. in both of these, x axis is time. rows should display activity of a set of channels. in the raster, color = magnitude, y axis = channel, x axis = time. default to ‘parula’ colormap; include a dropdown to choose from 15 different popular colormaps. in the ‘timecourse’ version, each row shows the timecourse of a channel as black on white background; within that row, y axis shows magnitude of that channel at each timepoint. 
% 
% use the ‘labels’ field of the fieldtrip object to label the rows corresponding to each channel. use ‘time’ field from the fieldtrip variable to label x axis. 
% 
% include parameter op.n_chans_raster (default to 40) and op.n_chans_timecourse (default to 20) which determine how many channels to include on the screen at once. 
% 
% The user should be able to select groups of [time x channel] rectangles. in raster mode, they should then be able to click and drag a portion of the raster, to select rectangular groups of coordinates (channel x timepoints). these should stay visibly selected with a gray border after the user has made a selection, so that the user can highlight multiple groups of coordinates. this should work basically the same for timecourses, except that portions of timecourse plots across 1 or more channels will get highlighted at once. When a group is selected in raster mode, it should be surrounded by a gray box. When selected in timecourse mode, the times x channels should be highlighted in light gray. 
% 
% when the gui is first opened, create a table variable called ‘artifact’ with the following columns:
% -id - double -  just a number indicating the row of the table [1 to nrows]
% -starts - double - time of the start of the artifactual window in global time coordinates [GTC] - seconds since midnight on the day of the experiment
% -ends - double - end of the artifact window, also in GTC
% -duration - double - starts minus ends
% -label - string - channel name - must correspond to a label in the ‘label’ field in the fieldtrip object
% 
% The ultimate output of this GUI is intended to be a .tsv table in this format. 
% 
% include a button ‘add selected artifact’. when this is clicked, then all grey selections should turn red and stay red unless they are removed. additionally, add all selected [times x channels] to the artifact table. This should add one row per channel selected; the starts, ends, and duration values in this row should indicate the time selected. however, if any of these [time x channels] are contiguous with a row already in the table, combine the two into a single row in the artifact table. once the rows have been added to the artifact table, sort the table - first by ‘starts’ then by ‘label’ to break ties. do this combination of artifact tables and updating of the plot with a subfunction ‘combine_artifact_tables’, which we will also use later. 
% 
% include a button for ‘remove selected artifact’. The user should be able to highlight [time x channel]s that have already been added to the artifact table (as part of the grey rectangle they are selecting). When ‘remove selected artifact’ is clicked, then all selected [time x channels] should be removed from the artifact table, and the red outline/highlight around them should be removed. if this would remove more than one non-contiguous groups of [time x channels], bring up a popup dialogue box asking “Remove multiple groups of time x channels?” with ‘yes’ and ‘no’ options. Do the removal if ‘yes’ is selected. If ‘no’ is selected, then don’t do the removal, but keep the gray selection box/highlight where it is. 
% 
% in the gui, include a dropdown, which lists chunks of trials [chunks sized at op.n_chans_raster or op.n_chans_timecourse]. when the user selects a chunk of trials, switch the display that chunk of channels. in the listing of channels in dropdown, list the first and last channel (e.g. “ecog_101 - dbs_52”). 
% 
% when switching between raster and timecourse mode, make sure to update the red outlines/highlights based on the current contents of the artifact table. also keep the gray ‘curent selection’ box where it is on the same [time x channel]s. The currently-viewed channels will likely not be the same, due to the different number of channels per view, so pick the channel group that overlaps the most with the previous view. some selected channels in the previous view may not be visible anymore, so remove those from the selected [time x channel]s. 
% 
% the gui should have an option ‘load artifact mask’. when this is pressed, if the artifact table is currently non-empty - if any [time x channel]s are highlighted in red - display a dialogue box saying “Artifact table is not empty. Add loaded table to current selection?” with ‘yes’ and ‘no’ options. this ‘load artifact mask’ option should open a dialogue to load a file. start by trying to look in the following folder: 
% Y:\DBS\derivatives\sub-[SUBJECT]\annot\
% 
% This should look for files with ‘artifact’ in the name which end with ‘.tsv’. First load this into matabl as a variable called ‘artifact_loaded’. It should have the same columns as described above in the ‘artifact’ table. Then, give an error dialogue box if any of the following are true (and display the problem) [do not output an actual matlab error, just display the dialogue box and delete the ‘artifact_loaded’ variable]:
% -table has zero rows
% -missing any of the 4 column variables, or they are an unexpected class
% -durations are not equal to starts minus ends
% -durations are zero or less or are nans
% -any of the ‘label’s do not match a label in the fieldtrip object
% -any of the time windows exist outside of the time windows indicated by the fieldtrip ‘time’ field [meaning they are outside the time range of what we have loaded from the fieldtrip file]
% 
% If there are additional columns, remove them from the table. If the table passes the checks, then add its columns to the current ‘artifact’ variable, via the ‘combine_artifact_tables’ subfunction described above, which should also update the plot. if no artifacts have been added in the current work session [nothing has been highlighted in red], this should function the same, because the loaded artifact table will be added to the zero-row artifact table. then delete the ‘artifact_loaded’ variable, because its info should be copied to the ‘artifact’ variable. 
% 
% Include a button ‘save artifact mask’. When this is clicked, open a file browser dialogue box. start in the same folder we load artifact tables from. the default file savename is:
% sub-[SUBJECT]_ses-intraop_task-smsl_artifact-manual.tsv
% 
% this save subfunction should always warn about overwriting - though this should already be included if it’s a normal windows file-saving dialogue box. 
% 
% put all buttons, fields, and menus stacked on the left side; rasters/timecourses to the right.


% output a modified version of this script [not word doc] with these changes:
% -disable to ability to select multiple boxes at once. once a selection is made, erase the previous [gray] selections.
% -allow a keyboard shortcut - when user presses “a”, the selected grey box becomes an artifact box. also change the label on the ‘add artifact’ button to “Add artifact (a)” to indicate this.
% -add keyboard shortcut “r” for remove selected artifact; update the button to add “ (r)" to the end of the button label.
% -always start in ‘timecourse’ view mode, not raster
% -change the expected end of the fieldtip file string from '_ses-intraop_task-smsl_ft-raw_trial.mat' to 'ses-intraop_task-smsl_ft-raw.mat'
% -add a button ‘zoom selection’ and ‘zoom full’. when ‘zoom selection’ is clicked, change the time we are looking at and the channels we are looking at to only the selected [time x channel]s.. ‘zoom full’ should zoom
% ---include keyboard shortcut 'q' for zoom selection, and shortcut 'w' for zoom full. add these shortcut keys to the button labels.
% 
% -when trying to load load an artifact table, add a check: look for a string at the beginning of the artifact table name which will be ‘sub-[SUBJECT]’. if [SUBJECT] isn’t the same as what we currently have loaded for op.sub, show both of them to the user in a dialogue box, and give option to proceed or to not load artifact mask.
% -add 2 fillable fields on the side of the gui that display the values of  op.n_chans_timecourse and  op.n_chans_raster. the user should be able to edit these values, then click an ‘update’ button next to them which changes these values. this should immediately change the number of channels currently being viewed, and also update what is shown in the dropdown menus to selecting chunks of channels.
% -change how the values of the raster are mapped to colors: they should always be relative to values within the channel itself, so that extreme values in one channel don’t make values in other channels all appear to be zero

% make the following changes:
% -add up and down buttons next to the dropdown for channel group. these should shift to the next/previous channel group, as if they were selected from the dropdown
% -replace the dropdown for view mode with a single button: Timecourse // Raster. This should switch view modes when clicked. Very clearly highlight the current mode by highlighting that side in green and bolding that text. include keyboard shortcut label (m)
% -change these shortcuts and update button labels:
% --- switch raster/timecourse = m
% --- load artifact mask = o
% --- save artifact mask = k