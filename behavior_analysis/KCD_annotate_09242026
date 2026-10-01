function annotate_disfluencies_video(initialVideo)
% ANNOTATE_DISFLUENCIES_VIDEO  GUI for scoring stuttering / disfluency events from video.
%
%   annotate_disfluencies_video              opens the GUI (load a video with the button)
%   annotate_disfluencies_video(videoPath)   opens the GUI and loads the given video
%
% LAYOUT
%   Left  : video | spectrogram | waveform | transcript strip | notes strip | overview
%   Right : Events list (play / edit text / change type / delete, playback speed)
%           and a Quick tips box
%   Bottom: navigation row (go to time, zoom, saved time frames), file row, status bar
%
% WORKFLOW
%   1. Load video.
%   2. Drag on the spectrogram or waveform to select a stretch of speech.
%   3. Press R / B / P to label it (repetition / block / prolongation).
%   4. Type the transcript / notes in the white boxes on the strips (Enter saves).
%   5. Drag the black handles on either edge of a selected block to resize it.
%   6. Space = play / pause (pause leaves the cursor where it stopped).
%   7. Ctrl+S to save. "Load annotations" can merge into or replace what you have.
%
% All annotation blocks and their text are drawn in black; the outline style
% (solid / dashed / dotted) and the label tell the event types apart.
%
% All event times come from the AUDIO sample clock (resolution = audio, not video
% frames). During playback the video frame is slaved to the audio position.
%
% Requires base MATLAB only (no toolboxes). Audio is read with AUDIOREAD; on some
% Linux setups you may need to load a container whose audio AUDIOREAD can decode.
%
% Every button / key / mouse callback is wrapped so an error pops up with the real
% message and the line number in this file.

% ----------------------------------------------------------------------------
% Configuration
% ----------------------------------------------------------------------------
APPNAME = 'Disfluency Annotator';
APPVER  = '2.0';

% Event types: name, trigger key, outline style, colour. Extend freely.
% 'plain' = a transcript / notes segment with no disfluency event: it gets no
% coloured block, only a thin grey outline on the Transcript and Notes strips.
eventTypes = struct( ...
    'name',  {'repetition', 'block', 'prolongation', 'text-only'}, ...
    'key',   {'r',          'b',     'p',            't'}, ...
    'style', {'-',          '--',    ':',            '-'}, ...
    'color', {[0.85 0.15 0.15], [0.15 0.35 0.90], [0.10 0.60 0.20], [0.55 0.55 0.55]}, ...
    'plain', {false,        false,   false,          true});
TXT_TYPE = find([eventTypes.plain],1);

ANN_COLOR  = [0 0 0];            % all annotation TEXT (labels, transcript, notes)
SEL_COLOR  = [0.45 0.45 0.45];   % active selection (neutral so type colours stay clear)
OV_COLOR   = [0.20 0.45 0.95];   % overview 'you are here' window
CUR_COLOR  = [0.90 0.50 0.00];   % cursor
PLAY_COLOR = [0.85 0.00 0.85];   % playhead
BG         = [0.94 0.94 0.94];
SPEEDS     = [1 0.75 0.5 0.25];

MINVIEW = 0.02;     % narrowest view (s)
MINSEL  = 0.005;    % shortest selection / event (s)
EDGE_PX = 8;        % grab tolerance for the resize handles (pixels)
MAXUNDO = 40;       % undo depth

% ----------------------------------------------------------------------------
% Shared state
% ----------------------------------------------------------------------------
initPath = '';
if nargin >= 1 && ~isempty(initialVideo), initPath = initialVideo; end

vid = []; videoPath = ''; audio = []; fs = 0; dur = 0; audioMax = 1; maxFreq = 5000;
viewStart = 0; viewEnd = 1; cursorTime = 0; selStart = NaN; selEnd = NaN;
events = struct('start',{},'end',{},'type',{},'transcript',{},'notes',{});
frames = struct('name',{},'a',{},'b',{});
currentEvent = []; dirty = false; undoStack = {}; listMap = [];

isEditing = false; editIdx = []; editPanel = ''; typingBox = []; transBox = []; notesBox = [];
isPlaying = false; playStartTime = 0; playEndTime = 0; playSpeed = 1; mainPlayer = [];
dragAxes = []; dragMode = ''; dragStartT = 0; downPix = [0 0]; didDrag = false;
dragPanel = ''; moveIdx = []; moveOrig = [0 0]; movePushed = false;
resizeUndoPushed = false;

imgVideo = []; placeholderTxt = []; hSpecImg = []; hWave = [];
evGfx = gobjects(0); selGfx = gobjects(0); curGfx = gobjects(0); playGfx = gobjects(0);
ovEnvGfx = gobjects(0); ovDyn = gobjects(0); ovCur = gobjects(0);
ovMode = ''; ovGrab = 0; edtFrom = []; edtTo = []; btnGoTF = [];

fig = []; statusTxt = []; infoTxt = []; viewTxt = []; lstEvents = [];
popFrames = []; popSpeed = []; videoDependent = [];
splashFig = []; splashBar = []; splashMsg = [];

% ----------------------------------------------------------------------------
% Loading screen
% ----------------------------------------------------------------------------
makeSplash();
splashCleanup = onCleanup(@closeSplash); %#ok<NASGU>
splashStep(0.10,'Preparing workspace...');

% ----------------------------------------------------------------------------
% Main window
% ----------------------------------------------------------------------------
scrSz = get(0,'ScreenSize');
figW = min(1280, scrSz(3)-60); figH = min(940, scrSz(4)-90);
fig = figure('Name',APPNAME,'NumberTitle','off','Units','pixels', ...
    'Position',[max(20,(scrSz(3)-figW)/2) max(40,(scrSz(4)-figH)/2) figW figH], ...
    'Color',BG,'MenuBar','none','Toolbar','none','Visible','off', ...
    'WindowKeyPressFcn',@(~,e)cb(@()onKey(e)), ...
    'WindowScrollWheelFcn',@(~,e)cb(@()onScroll(e)), ...
    'WindowButtonMotionFcn',@(~,~)onHover(), ...
    'SizeChangedFcn',@(~,~)onResize(), ...
    'CloseRequestFcn',@(~,~)onClose());

splashStep(0.30,'Building video and signal panels...');

axVideo = axes('Parent',fig,'Units','normalized','Position',[0.04 0.60 0.70 0.37], ...
    'Color',[0 0 0]);
axis(axVideo,[0 1 0 1]); axis(axVideo,'off');
placeholderTxt = text(axVideo,0.5,0.5,'Click "Load video" to begin', ...
    'HorizontalAlignment','center','FontSize',14,'Color',[0.6 0.6 0.6]);

axSpec     = axes('Parent',fig,'Units','normalized','Position',[0.07 0.440 0.67 0.145]);
axWave     = axes('Parent',fig,'Units','normalized','Position',[0.07 0.355 0.67 0.075]);
axTrans    = axes('Parent',fig,'Units','normalized','Position',[0.07 0.290 0.67 0.055]);
axNotes    = axes('Parent',fig,'Units','normalized','Position',[0.07 0.225 0.67 0.055]);
axOverview = axes('Parent',fig,'Units','normalized','Position',[0.07 0.138 0.67 0.036]);

for axInit = [axSpec axWave axTrans axNotes axOverview]
    hold(axInit,'on'); box(axInit,'on');
    set(axInit,'XTick',[],'YTick',[],'Color',[1 1 1], ...
        'XColor',[0.15 0.15 0.15],'YColor',[0.15 0.15 0.15]);
end
set(axSpec,'YDir','normal','YTickMode','auto'); ylabel(axSpec,'Freq (Hz)');
ylabel(axWave,'Amp'); ylim(axWave,[-1 1]);
ylabel(axTrans,'Transcript','Rotation',0,'HorizontalAlignment','right'); ylim(axTrans,[0 1]);
ylabel(axNotes,'Notes','Rotation',0,'HorizontalAlignment','right');      ylim(axNotes,[0 1]);
ylabel(axOverview,'Overview','Rotation',0,'HorizontalAlignment','right'); ylim(axOverview,[0 1]);

set(axSpec,    'ButtonDownFcn',@(~,~)cb(@()onAxDown(axSpec)));
set(axWave,    'ButtonDownFcn',@(~,~)cb(@()onAxDown(axWave)));
set(axTrans,   'ButtonDownFcn',@(~,~)cb(@()onAxDown(axTrans,'transcript')));
set(axNotes,   'ButtonDownFcn',@(~,~)cb(@()onAxDown(axNotes,'notes')));
set(axOverview,'ButtonDownFcn',@(~,~)cb(@onOverviewDown));

% Inline text boxes shown over the selected event on the Transcript / Notes
% strips. Click one to type; Enter (or clicking anywhere else) saves it.
transBox = uicontrol(fig,'Style','edit','Max',1,'Min',0,'Units','pixels', ...
    'HorizontalAlignment','left','FontSize',10,'Visible','off','Enable','inactive', ...
    'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Tag','transcript', ...
    'TooltipString','Click to type the transcript - Enter saves', ...
    'ButtonDownFcn',@(~,~)cb(@()beginTyping('transcript')), ...
    'Callback',@(src,~)cb(@()commitBox(src)));
notesBox = uicontrol(fig,'Style','edit','Max',1,'Min',0,'Units','pixels', ...
    'HorizontalAlignment','left','FontSize',10,'Visible','off','Enable','inactive', ...
    'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Tag','notes', ...
    'TooltipString','Click to type notes - Enter saves', ...
    'ButtonDownFcn',@(~,~)cb(@()beginTyping('notes')), ...
    'Callback',@(src,~)cb(@()commitBox(src)));


splashStep(0.55,'Creating controls...');

% ---- right column: events panel ----------------------------------------------
pnlEv = uipanel(fig,'Title','Events','Units','normalized', ...
    'Position',[0.765 0.30 0.22 0.675],'BackgroundColor',BG,'FontWeight','bold');
% Colour-coded events table: each row is tinted with its event colour;
% text-only segments stay white.
lstEvents = uitable(pnlEv,'Units','normalized','Position',[0.04 0.40 0.92 0.58], ...
    'ColumnName',{'#','Start','Len (s)','Type','Transcript'}, ...
    'ColumnWidth',{30 70 50 46 110},'ColumnEditable',false(1,5),'RowName',[], ...
    'Data',cell(0,5),'FontSize',9,'ForegroundColor',[0 0 0], ...
    'BackgroundColor',[1 1 1],'RowStriping','on', ...
    'CellSelectionCallback',@(src,e)cb(@()onListSelect(src,e)));
btnPlayEv = mkBtn(pnlEv,'Play event',  [0.04 0.320 0.45 0.065],@playCurrentEvent, ...
    'Play the selected event');
btnEditEv = mkBtn(pnlEv,'Edit text',   [0.51 0.320 0.45 0.065],@editCurrentEvent, ...
    'Type the transcript of the selected event');
btnTypeEv = mkBtn(pnlEv,'Change type', [0.04 0.245 0.45 0.065],@changeTypeDialog, ...
    'Change the event type of the selected event');
btnDelEv  = mkBtn(pnlEv,'Delete event',[0.51 0.245 0.45 0.065],@deleteCurrentEvent, ...
    'Delete the selected event (Ctrl+Z to undo)');
btnTextEv = mkBtn(pnlEv,'Add transcript / notes only  (T)',[0.04 0.170 0.92 0.065], ...
    @()applyType(TXT_TYPE),'Add transcript / notes for the selection without making it an event');
uicontrol(pnlEv,'Style','text','String','Playback speed:','Units','normalized', ...
    'Position',[0.04 0.100 0.45 0.045],'HorizontalAlignment','left','BackgroundColor',BG);
popSpeed = uicontrol(pnlEv,'Style','popupmenu','Units','normalized', ...
    'String',{'1.00x (normal)','0.75x','0.50x','0.25x'},'Value',1, ...
    'Position',[0.51 0.105 0.45 0.050],'Callback',@(src,~)cb(@()onSpeed(src)));
infoTxt = uicontrol(pnlEv,'Style','text','String','No selection','Units','normalized', ...
    'Position',[0.04 0.010 0.92 0.085],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontSize',8.5);

pnlTips = uipanel(fig,'Title','Quick tips','Units','normalized', ...
    'Position',[0.765 0.135 0.22 0.155],'BackgroundColor',BG,'FontWeight','bold');
uicontrol(pnlTips,'Style','text','Units','normalized','Position',[0.03 0.02 0.94 0.96], ...
    'HorizontalAlignment','left','BackgroundColor',BG,'FontSize',8.5,'String', { ...
    '1. Drag on the spectrogram / waveform to select.', ...
    '2. Press R, B or P to label it - or T for text only.', ...
    '3. Drag the black edge handles to resize.', ...
    '4. Space = play / pause.', ...
    '5. Click a block, then a strip, to edit text.', ...
    'F1 or Help = full instructions'});

% ---- navigation row ------------------------------------------------------------
uicontrol(fig,'Style','text','String','Navigate:','Units','normalized', ...
    'Position',[0.04 0.083 0.055 0.035],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontWeight','bold');
btnGoto  = mkBtn(fig,'Go to time...',  [0.095 0.087 0.085 0.040],@goToTimeDialog, ...
    'Type a start and end time to look at (Ctrl+G)');
btnZin   = mkBtn(fig,'Zoom in',        [0.185 0.087 0.065 0.040],@()zoomAbout(cursorCenter(),0.5), ...
    'Zoom in around the cursor (Ctrl+I)');
btnZout  = mkBtn(fig,'Zoom out',       [0.255 0.087 0.065 0.040],@()zoomAbout(cursorCenter(),2), ...
    'Zoom out around the cursor (Ctrl+O)');
btnZsel  = mkBtn(fig,'Fit selection',  [0.325 0.087 0.080 0.040],@zoomToSelection, ...
    'Zoom to the current selection (Ctrl+N)');
btnZall  = mkBtn(fig,'Full recording', [0.410 0.087 0.085 0.040],@fullView, ...
    'Show the whole recording (Ctrl+A)');
popFrames = uicontrol(fig,'Style','popupmenu','String',{'Saved time frames...'}, ...
    'Units','normalized','Position',[0.500 0.087 0.170 0.040], ...
    'TooltipString','Jump to a saved time frame', ...
    'Callback',@(src,~)cb(@()onFramePick(src)));
btnFsave = mkBtn(fig,'Save frame',     [0.675 0.087 0.075 0.040],@saveFrame, ...
    'Save the current view as a named time frame');
btnFdel  = mkBtn(fig,'Remove frame',   [0.755 0.087 0.080 0.040],@removeFrame, ...
    'Remove the time frame chosen in the drop-down');
% current-view readout, sitting just above the Overview strip
viewTxt = uicontrol(fig,'Style','text','String','','Units','normalized', ...
    'Position',[0.07 0.176 0.36 0.024],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontSize',9);
% type an exact time frame (seconds or m:ss); Enter or Go applies it
uicontrol(fig,'Style','text','String','Show from','Units','normalized', ...
    'Position',[0.432 0.176 0.055 0.024],'HorizontalAlignment','right', ...
    'BackgroundColor',BG,'FontSize',9);
edtFrom = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.490 0.174 0.085 0.028], ...
    'FontSize',9,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Enable','inactive', ...
    'Tag','tfrom','HorizontalAlignment','center', ...
    'TooltipString','Start time (seconds or m:ss) - click, type, Enter', ...
    'ButtonDownFcn',@(src,~)cb(@()beginFieldEdit(src)),'Callback',@(src,~)cb(@()commitBox(src)));
uicontrol(fig,'Style','text','String','to','Units','normalized', ...
    'Position',[0.577 0.176 0.018 0.024],'HorizontalAlignment','center', ...
    'BackgroundColor',BG,'FontSize',9);
edtTo = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.597 0.174 0.085 0.028], ...
    'FontSize',9,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Enable','inactive', ...
    'Tag','tto','HorizontalAlignment','center', ...
    'TooltipString','End time (seconds or m:ss) - click, type, Enter', ...
    'ButtonDownFcn',@(src,~)cb(@()beginFieldEdit(src)),'Callback',@(src,~)cb(@()commitBox(src)));
btnGoTF = mkBtn(fig,'Go',[0.686 0.173 0.054 0.030],@applyTimeFields, ...
    'Show the time frame typed in the From / To boxes');

% ---- file row --------------------------------------------------------------------
mkBtn(fig,'Load video',            [0.040 0.035 0.100 0.042],@uiLoadVideo, ...
    'Open a video file');
btnLoadAnnot = mkBtn(fig,'Load annotations...',[0.145 0.035 0.120 0.042],@loadAnnotations, ...
    'Merge an annotation file into the current one, or replace it');
btnSave  = mkBtn(fig,'Save annotations',[0.270 0.035 0.120 0.042],@saveAnnotations, ...
    'Save all events to an Excel file (Ctrl+S)');
btnUndo  = mkBtn(fig,'Undo',            [0.395 0.035 0.070 0.042],@undo, ...
    'Undo the last change (Ctrl+Z)');
mkBtn(fig,'Help',                  [0.470 0.035 0.070 0.042],@showHelp, ...
    'Full instructions (F1)');
mkBtn(fig,'Shortcuts',             [0.545 0.035 0.080 0.042],@showShortcuts, ...
    'Keyboard and mouse shortcuts');

% colour legend: swatch + black label per event type
legX = 0.635;
for tiInit = 1:numel(eventTypes)
    swatch = eventTypes(tiInit).color;
    if eventTypes(tiInit).plain, swatch = [0.88 0.88 0.88]; end
    uicontrol(fig,'Style','text','String','','Units','normalized', ...
        'Position',[legX 0.047 0.012 0.018],'BackgroundColor',swatch);
    uicontrol(fig,'Style','text','Units','normalized','Position',[legX+0.015 0.038 0.074 0.032], ...
        'String',sprintf('%s = %s',upper(eventTypes(tiInit).key),eventTypes(tiInit).name), ...
        'HorizontalAlignment','left','BackgroundColor',BG,'FontSize',8.5);
    legX = legX + 0.089;
end

statusTxt = uicontrol(fig,'Style','text','String','No video loaded - click "Load video" to begin.', ...
    'Units','normalized','Position',[0.04 0.004 0.94 0.026],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontSize',9,'ForegroundColor',[0.2 0.2 0.2]);

% Force black label text and panel titles (macOS dark mode otherwise draws
% them in a pale grey that is hard to read on the light background).
set(findall(fig,'Style','text'),'ForegroundColor',[0 0 0]);
set(findall(fig,'Type','uipanel'),'ForegroundColor',[0 0 0]);

videoDependent = [edtFrom edtTo btnGoTF lstEvents btnPlayEv btnEditEv btnTypeEv btnDelEv btnTextEv popSpeed ...
    btnGoto btnZin btnZout btnZsel btnZall popFrames btnFsave btnFdel ...
    btnLoadAnnot btnSave btnUndo];
set(videoDependent,'Enable','off');

splashStep(0.80,'Setting up playback engine...');
playTimer = timer('ExecutionMode','fixedRate','Period',0.03,'BusyMode','drop', ...
    'TimerFcn',@(~,~)cb(@playTick));

splashStep(1.00,'Ready');
closeSplash();
set(fig,'Visible','on'); drawnow;

showWelcome();
if ~isempty(initPath), cb(@()loadVideo(initPath)); end

% ============================================================================
%                               NESTED FUNCTIONS
% ============================================================================

    % ---------------- universal callback guard ------------------------------
    function cb(fn)
        try
            fn();
        catch err
            loc = '';
            frame = [];
            for si = 1:numel(err.stack)
                if strcmp(err.stack(si).name,'annotate_disfluencies_video') || ...
                        contains(err.stack(si).name,'annotate_disfluencies_video/')
                    frame = err.stack(si); break;
                end
            end
            if isempty(frame) && ~isempty(err.stack), frame = err.stack(1); end
            if ~isempty(frame)
                loc = sprintf('\n\n(in %s, line %d)',frame.name,frame.line);
            end
            if ~isempty(statusTxt) && ishghandle(statusTxt)
                set(statusTxt,'String',['ERROR: ' err.message],'ForegroundColor',[0.8 0 0]);
            end
            errordlg([err.message loc],'Error');
        end
    end

    % ---------------- buttons -----------------------------------------------
    function h = mkBtn(parent,label,pos,fn,tip)
        h = uicontrol(parent,'Style','pushbutton','String',label,'Units','normalized', ...
            'Position',pos,'TooltipString',tip, ...
            'Callback',@(src,~)cb(@()btnPress(src,fn)));
    end

    function btnPress(src,fn)
        % Toggling Enable off/on removes keyboard focus from the button, so
        % keys (Space, R/B/P ...) go to the figure and do not re-click it.
        try
            set(src,'Enable','off'); drawnow; set(src,'Enable','on');
        catch
        end
        if isEditing, commitEdit(); end
        fn();
    end

    % ---------------- loading screen ----------------------------------------
    function makeSplash()
        ss = get(0,'ScreenSize'); sw = 460; sh = 220;
        splashFig = figure('Name',[APPNAME ' - loading'],'NumberTitle','off', ...
            'MenuBar','none','ToolBar','none','Resize','off','Color',[1 1 1], ...
            'Units','pixels','Position',[(ss(3)-sw)/2 (ss(4)-sh)/2 sw sh], ...
            'HandleVisibility','off','CloseRequestFcn','');
        sax = axes('Parent',splashFig,'Units','normalized','Position',[0 0 1 1], ...
            'XLim',[0 1],'YLim',[0 1],'Visible','off','HandleVisibility','off');
        text(sax,0.5,0.80,APPNAME,'FontSize',18,'FontWeight','bold', ...
            'HorizontalAlignment','center');
        text(sax,0.5,0.66,['Scoring stuttering events from video  |  v' APPVER], ...
            'FontSize',10,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');
        patch(sax,[0.1 0.9 0.9 0.1],[0.40 0.40 0.47 0.47],[0.92 0.92 0.92], ...
            'EdgeColor',[0.6 0.6 0.6]);
        splashBar = patch(sax,[0.1 0.1 0.1 0.1],[0.40 0.40 0.47 0.47],[0.1 0.1 0.1], ...
            'EdgeColor','none');
        splashMsg = text(sax,0.5,0.27,'Starting...','FontSize',9, ...
            'HorizontalAlignment','center','Color',[0.25 0.25 0.25]);
        drawnow;
    end

    function splashStep(frac,msg)
        if isempty(splashFig) || ~ishghandle(splashFig), return; end
        x1 = 0.1 + 0.8*frac;
        set(splashBar,'XData',[0.1 x1 x1 0.1]); set(splashMsg,'String',msg);
        drawnow; pause(0.15);
    end

    function closeSplash()
        if ~isempty(splashFig) && ishghandle(splashFig), delete(splashFig); end
        splashFig = [];
    end

    function showWelcome()
        try
            if ispref('DisfluencyAnnotator','hideWelcome') && ...
                    getpref('DisfluencyAnnotator','hideWelcome')
                return;
            end
        catch
        end
        msg = sprintf([ ...
            'Quick start:\n\n' ...
            '1. Click "Load video" and choose a recording.\n' ...
            '2. Click and drag on the spectrogram or waveform to select a stretch of speech.\n' ...
            '3. Press R (repetition), B (block) or P (prolongation) to label it.\n' ...
            '4. Type the transcript and notes in the white boxes on the strips (Enter saves).\n' ...
            '5. Drag the black handles on either edge of a block to fine-tune it.\n' ...
            '6. Press Space to play and Space again to pause.\n' ...
            '7. Save often with Ctrl+S.\n\n' ...
            'Press F1 or the Help button at any time for full instructions.']);
        c = questdlg(msg,['Welcome to ' APPNAME],'Get started','Don''t show again','Get started');
        if strcmp(c,'Don''t show again')
            try, setpref('DisfluencyAnnotator','hideWelcome',true); catch, end
        end
    end

    % ---------------- loading video -----------------------------------------
    function uiLoadVideo()
        if ~confirmUnsaved('loading another video'), return; end
        [fn,fp] = uigetfile( ...
            {'*.mp4;*.mov;*.avi;*.m4v;*.mkv;*.wmv','Video files';'*.*','All files'}, ...
            'Select a video file');
        if isequal(fn,0), return; end
        loadVideo(fullfile(fp,fn));
    end

    function loadVideo(path)
        wb = waitbar(0,'Opening video file...','Name','Loading video');
        wbClean = onCleanup(@()deleteValid(wb)); %#ok<NASGU>
        try
            v = VideoReader(path);
        catch err
            deleteValid(wb);
            errordlg(sprintf('Could not open the video:\n\n%s',err.message),'Load error'); return;
        end
        waitbar(0.25,wb,'Extracting the audio track (long videos take a moment)...');
        try
            [a,fsr] = audioread(path);
        catch
            deleteValid(wb);
            errordlg(sprintf(['Could not read the audio from this file.\n\nOn some systems ' ...
                'AUDIOREAD cannot decode compressed video audio. Convert the video ' ...
                '(e.g. to .mp4 with AAC audio) or extract a .wav and try again.']), ...
                'Audio error'); return;
        end
        if isempty(a)
            deleteValid(wb);
            errordlg('This video has no audio track, so it cannot be annotated here.','Audio error');
            return;
        end
        waitbar(0.55,wb,'Preparing audio...');
        stopPlay(); commitEditQuiet();
        vid = v; videoPath = path; fs = fsr;
        if size(a,2) > 1, a = mean(a,2); end
        audio = single(a(:)); audioMax = max(max(abs(audio)),eps);
        dur = numel(audio)/fs; maxFreq = min(5000, fs/2);
        events = struct('start',{},'end',{},'type',{},'transcript',{},'notes',{});
        frames = struct('name',{},'a',{},'b',{});
        currentEvent = []; selStart = NaN; selEnd = NaN; cursorTime = 0;
        undoStack = {}; dirty = false;
        viewStart = 0; viewEnd = min(dur,10);
        buildPlayer();

        waitbar(0.70,wb,'Building the recording overview...');
        buildOverview();

        waitbar(0.85,wb,'Rendering the first frame...');
        if ~isempty(placeholderTxt) && ishghandle(placeholderTxt)
            delete(placeholderTxt); placeholderTxt = [];
        end
        cla(axVideo);
        try
            vid.CurrentTime = 0; frame0 = readFrame(vid);
        catch
            frame0 = zeros(2,2,3,'uint8');
        end
        imgVideo = image(axVideo,frame0); axis(axVideo,'image','off');
        set(imgVideo,'PickableParts','none','HitTest','off');

        set(videoDependent,'Enable','on');
        updateFramesPopup(0); updateEventList(); updateTitle();
        waitbar(1,wb,'Done');
        refreshView(); showFrameAt(0);
        [~,nm,ext] = fileparts(path);
        setStatus(sprintf(['Loaded %s%s  (%s, %d Hz audio, %.0f fps).  Drag on the ' ...
            'spectrogram to select speech, then press R / B / P.'], ...
            nm,ext,fmtTime(dur),fs,vid.FrameRate),'ok');
    end

    % ---------------- view / drawing ----------------------------------------
    function refreshView()
        if isempty(vid), return; end
        viewStart = max(0,viewStart); viewEnd = min(dur,viewEnd);
        if viewEnd - viewStart < MINVIEW, viewEnd = min(dur,viewStart+MINVIEW); end
        set([axSpec axWave axTrans axNotes],'XLim',[viewStart viewEnd]);
        set([axSpec axWave axTrans],'XTick',[]); set(axNotes,'XTickMode','auto');
        safeDraw(@drawSpectrogram,'spectrogram');
        safeDraw(@drawWaveform,'waveform');
        safeDraw(@redrawEvents,'events');
        updateSelectionGraphics(); updateCursorGraphics();
        set(viewTxt,'String',sprintf('View: %s - %s  (%.2f s wide)', ...
            fmtTime(viewStart),fmtTime(viewEnd),viewEnd-viewStart));
        syncTimeFields(viewStart,viewEnd);
    end

    function safeDraw(fn,label)
        try
            fn();
        catch err
            setStatus(['Draw error (' label '): ' err.message],'warn');
        end
    end

    function drawSpectrogram()
        i0 = max(1,floor(viewStart*fs)+1); i1 = min(numel(audio),ceil(viewEnd*fs));
        if ~isempty(hSpecImg)&&ishghandle(hSpecImg), delete(hSpecImg); hSpecImg=[]; end
        if i1-i0 < 64, return; end
        seg = double(audio(i0:i1));
        winlen = max(64,round(0.006*fs)); if mod(winlen,2)==0, winlen=winlen+1; end
        Ncol = 700;
        hop  = max(1, floor((numel(seg)-winlen)/Ncol));
        nfft = 2^nextpow2(max(winlen,512));
        w    = 0.5 - 0.5*cos(2*pi*(0:winlen-1)'/(winlen-1));   % Hann window (no toolbox)
        st   = 1:hop:(numel(seg)-winlen+1);
        if isempty(st), return; end
        nf = floor(nfft/2)+1;
        P  = zeros(nf,numel(st));
        for c = 1:numel(st)
            idx = st(c):st(c)+winlen-1;
            X = fft(seg(idx).*w, nfft);
            P(:,c) = abs(X(1:nf));
        end
        P = 20*log10(P+eps);
        F = (0:nf-1)'*(fs/nfft);
        Tc = (st-1+winlen/2)/fs + viewStart;
        fmask = F <= maxFreq; if ~any(fmask), fmask = true(size(F)); end
        hSpecImg = imagesc(axSpec, Tc, F(fmask), P(fmask,:));
        set(axSpec,'YDir','normal','YLim',[0 maxFreq],'XLim',[viewStart viewEnd]);
        mx = max(P(fmask,:),[],'all');
        if isfinite(mx), setCLim(axSpec,[mx-70 mx]); end
        colormap(axSpec, flipud(gray(256)));
        set(hSpecImg,'PickableParts','none','HitTest','off'); uistack(hSpecImg,'bottom');
    end

    function setCLim(ax,lims)
        try, clim(ax,lims); catch, caxis(ax,lims); end %#ok<CAXIS>
    end

    function drawWaveform()
        i0 = max(1,floor(viewStart*fs)+1); i1 = min(numel(audio),ceil(viewEnd*fs));
        if ~isempty(hWave), delete(hWave(ishghandle(hWave))); hWave = []; end
        if i1 <= i0, return; end
        seg = double(audio(i0:i1)); tt = ((i0:i1)-1)/fs;
        if numel(seg) <= 4000
            hWave = plot(axWave,tt,seg/audioMax,'Color',[0 0 0]);
        else
            nb = 1500; edges = round(linspace(1,numel(seg)+1,nb+1));
            mins = zeros(1,nb); maxs = zeros(1,nb); tc = zeros(1,nb);
            for k = 1:nb
                idx = edges(k):edges(k+1)-1; if isempty(idx), idx = edges(k); end
                sk = seg(idx); mins(k)=min(sk); maxs(k)=max(sk); tc(k)=tt(idx(1));
            end
            xp = [tc fliplr(tc)]; yp = [maxs fliplr(mins)]/audioMax;
            hWave = patch(axWave,'XData',xp,'YData',yp,'FaceColor',[0.1 0.1 0.1], ...
                'EdgeColor',[0.1 0.1 0.1]);
        end
        set(hWave,'PickableParts','none','HitTest','off');
        set(axWave,'XLim',[viewStart viewEnd],'YLim',[-1 1]); uistack(hWave,'bottom');
    end

    function redrawEvents()
        deleteValid(evGfx); evGfx = gobjects(0);
        vw = viewEnd - viewStart;
        yl = get(axSpec,'YLim');
        for k = 1:numel(events)
            s = events(k).start; en = events(k).end;
            if en < viewStart || s > viewEnd, continue; end
            ls  = eventTypes(events(k).type).style;
            col = eventTypes(events(k).type).color;
            isPlain = eventTypes(events(k).type).plain;
            isSel = ~isempty(currentEvent) && currentEvent == k;
            if isSel, fa = 0.32; lw = 2.5; else, fa = 0.18; lw = 1.5; end
            xs = [s en en s];
            tx0 = max(s,viewStart) + 0.004*vw;
            % stack text when a text-only segment overlaps a real event
            others = [1:k-1 k+1:numel(events)];
            mixed = false;
            for o = others
                if events(o).start < en && events(o).end > s && ...
                        eventTypes(events(o).type).plain ~= isPlain
                    mixed = true; break;
                end
            end
            if ~mixed, ty = 0.5; elseif isPlain, ty = 0.75; else, ty = 0.28; end
            if isPlain
                if isSel, ec = [0 0 0]; plw = 2; else, ec = col; plw = 1; end
                pTr = patch(axTrans,'XData',xs,'YData',[0.02 0.02 0.98 0.98], ...
                    'FaceColor','none','EdgeColor',ec,'LineWidth',plw);
                pNo = patch(axNotes,'XData',xs,'YData',[0.02 0.02 0.98 0.98], ...
                    'FaceColor','none','EdgeColor',ec,'LineWidth',plw);
                tTr = text(axTrans,tx0,ty,events(k).transcript,'Clipping','on', ...
                    'Interpreter','none','FontSize',9,'FontWeight','bold','Color',ANN_COLOR, ...
                    'VerticalAlignment','middle');
                tNo = text(axNotes,tx0,ty,events(k).notes,'Clipping','on', ...
                    'Interpreter','none','FontSize',9,'Color',ANN_COLOR, ...
                    'VerticalAlignment','middle');
                evGfx = [evGfx pTr pNo tTr tNo]; %#ok<AGROW>
                continue;
            end
            pSpec = patch(axSpec,'XData',xs,'YData',[yl(1) yl(1) yl(2) yl(2)], ...
                'FaceColor',col,'FaceAlpha',fa,'EdgeColor',col,'LineStyle',ls,'LineWidth',lw);
            pWave = patch(axWave,'XData',xs,'YData',[-1 -1 1 1], ...
                'FaceColor',col,'FaceAlpha',fa,'EdgeColor',col,'LineStyle',ls,'LineWidth',lw);
            pTr = patch(axTrans,'XData',xs,'YData',[0 0 1 1], ...
                'FaceColor',col,'FaceAlpha',0.15,'EdgeColor',col,'LineStyle',ls,'LineWidth',lw);
            pNo = patch(axNotes,'XData',xs,'YData',[0 0 1 1], ...
                'FaceColor',col,'FaceAlpha',0.15,'EdgeColor',col,'LineStyle',ls,'LineWidth',lw);
            pos = find(listMap == k,1);
            if isempty(pos), lab = upper(eventTypes(events(k).type).name);
            else, lab = sprintf('#%d %s',pos,upper(eventTypes(events(k).type).name)); end
            tLab = text(axSpec,tx0,yl(2)-0.03*diff(yl),lab,'Clipping','on', ...
                'Interpreter','none','FontSize',8,'FontWeight','bold','Color',ANN_COLOR, ...
                'BackgroundColor',[1 1 1],'Margin',1,'VerticalAlignment','top');
            tTr = text(axTrans,tx0,ty,events(k).transcript,'Clipping','on', ...
                'Interpreter','none','FontSize',9,'FontWeight','bold','Color',ANN_COLOR, ...
                'VerticalAlignment','middle');
            tNo = text(axNotes,tx0,ty,events(k).notes,'Clipping','on', ...
                'Interpreter','none','FontSize',9,'Color',ANN_COLOR, ...
                'VerticalAlignment','middle');
            evGfx = [evGfx pSpec pWave pTr pNo tLab tTr tNo]; %#ok<AGROW>
        end
        if ~isempty(evGfx), set(evGfx,'PickableParts','none','HitTest','off'); end
        drawOverview();
        positionBoxes();
    end

    function updateSelectionGraphics()
        deleteValid(selGfx); selGfx = gobjects(0);
        updateInfo();
        if isnan(selStart) || isnan(selEnd) || selEnd <= selStart || isempty(vid), return; end
        for hx = [axSpec axWave axTrans axNotes]
            yl = get(hx,'YLim');
            p = patch(hx,'XData',[selStart selEnd selEnd selStart], ...
                'YData',[yl(1) yl(1) yl(2) yl(2)],'FaceColor',SEL_COLOR, ...
                'FaceAlpha',0.12,'EdgeColor','none');
            selGfx = [selGfx p drawHandles(hx,selStart,selEnd)]; %#ok<AGROW>
        end
        set(selGfx,'PickableParts','none','HitTest','off'); stackEach(selGfx,'top');
    end

    function h = drawHandles(ax,t0,t1)
        % The same black edge "sliders" everywhere: a full-height line plus a
        % white grip, sized in pixels so they look identical on every strip.
        yl = get(ax,'YLim'); xl = get(ax,'XLim'); px = getpixelposition(ax,true);
        tpp = diff(xl)/max(px(3),1); ypp = diff(yl)/max(px(4),1);
        gw = 4*tpp; gh = min(0.6*diff(yl),24*ypp); ym = mean(yl);
        gy = [ym-gh/2 ym-gh/2 ym+gh/2 ym+gh/2];
        h = gobjects(0);
        for tt = [t0 t1]
            ln = plot(ax,[tt tt],yl,'-','Color',[0 0 0],'LineWidth',2.5);
            gp = patch(ax,'XData',tt+[-gw gw gw -gw],'YData',gy, ...
                'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
            h = [h ln gp]; %#ok<AGROW>
        end
        set(h,'PickableParts','none','HitTest','off');
    end

    function updateCursorGraphics()
        deleteValid(curGfx); curGfx = gobjects(0);
        if ~isempty(ovCur) && all(isgraphics(ovCur)), set(ovCur,'XData',[cursorTime cursorTime]); end
        if isnan(cursorTime) || isempty(vid), updateInfo(); return; end
        l1 = plot(axSpec,[cursorTime cursorTime],get(axSpec,'YLim'),'--','Color',CUR_COLOR,'LineWidth',1);
        l2 = plot(axWave,[cursorTime cursorTime],[-1 1],'--','Color',CUR_COLOR,'LineWidth',1);
        curGfx = [l1 l2];
        set(curGfx,'PickableParts','none','HitTest','off'); stackEach(curGfx,'top');
        updateInfo();
    end

    function updateInfo()
        if isempty(infoTxt) || ~ishghandle(infoTxt), return; end
        if isnan(selStart) || isnan(selEnd) || selEnd <= selStart
            s = {sprintf('Cursor: %s',fmtTime(cursorTime)),'No selection - drag on the timeline.'};
        else
            s = {sprintf('Selection: %s - %s',fmtTime(selStart),fmtTime(selEnd)), ...
                 sprintf('Length: %.3f s',selEnd-selStart)};
            if ~isempty(currentEvent) && currentEvent <= numel(events)
                s{2} = [s{2} sprintf('   |   %s',eventLabelShort(currentEvent))];
            end
        end
        set(infoTxt,'String',s);
    end

    % ---------------- overview strip ----------------------------------------
    function buildOverview()
        deleteValid(ovEnvGfx); ovEnvGfx = gobjects(0);
        n = numel(audio); nb = min(2000,n);
        L = floor(n/nb)*nb;
        M = reshape(abs(audio(1:L)),[],nb);
        env = double(max(M,[],1))/audioMax;
        tc = ((0:nb-1)+0.5)*(L/nb)/fs;
        set(axOverview,'XLim',[0 max(dur,eps)],'YLim',[0 1],'XTick',[]);
        ovEnvGfx = patch(axOverview,'XData',[tc fliplr(tc)], ...
            'YData',[0.5+0.45*env fliplr(0.5-0.45*env)], ...
            'FaceColor',[0.55 0.55 0.55],'EdgeColor','none');
        set(ovEnvGfx,'PickableParts','none','HitTest','off');
    end

    function drawOverview()
        deleteValid(ovDyn); ovDyn = gobjects(0); ovCur = gobjects(0);
        if isempty(vid), return; end
        pv = patch(axOverview,'XData',[viewStart viewEnd viewEnd viewStart],'YData',[0 0 1 1], ...
            'FaceColor',OV_COLOR,'FaceAlpha',0.25,'EdgeColor',OV_COLOR,'LineWidth',1.5);
        ovDyn = [ovDyn pv];
        for q = 1:numel(eventTypes)
            if isempty(events), break; end
            if eventTypes(q).plain, continue; end
            ev = events([events.type] == q);
            if isempty(ev), continue; end
            xs = reshape([[ev.start]; [ev.end]; nan(1,numel(ev))],1,[]);
            pe = plot(axOverview,xs,0.12*ones(size(xs)),'-','Color',eventTypes(q).color,'LineWidth',4);
            ovDyn = [ovDyn pe]; %#ok<AGROW>
        end
        % left / right handles ("sliders") on the view window
        ovDyn = [ovDyn drawHandles(axOverview,viewStart,viewEnd)];
        ovCur = plot(axOverview,[cursorTime cursorTime],[0 1],'-','Color',CUR_COLOR,'LineWidth',1.2);
        ovDyn = [ovDyn ovCur];
        set(ovDyn,'PickableParts','none','HitTest','off');
    end

    function onOverviewDown()
        if isempty(vid), return; end
        commitEdit();
        cp = get(axOverview,'CurrentPoint'); t = clampT(cp(1,1));
        ovMode = ovHit(t);
        if isempty(ovMode)                 % outside the window: centre on click
            w = viewEnd - viewStart; setView(t-w/2,t+w/2);
            ovMode = 'pan';
        end
        ovGrab = t - viewStart;
        set(fig,'WindowButtonMotionFcn',@(~,~)cb(@ovDrag), ...
            'WindowButtonUpFcn',@(~,~)cb(@endOverviewDrag));
    end

    function ovDrag()
        cp = get(axOverview,'CurrentPoint'); t = clampT(cp(1,1));
        switch ovMode
            case 'left',  setView(min(t,viewEnd-MINVIEW),viewEnd);
            case 'right', setView(viewStart,max(t,viewStart+MINVIEW));
            otherwise
                w = viewEnd - viewStart; a = t - ovGrab; setView(a,a+w);
        end
        drawnow limitrate;
    end

    function m = ovHit(t)
        % 'left' / 'right' near a handle, 'pan' inside the window, '' outside
        m = '';
        tol = 6*ovTimePerPixel();
        if abs(t-viewStart) <= tol, m = 'left';
        elseif abs(t-viewEnd) <= tol, m = 'right';
        elseif t > viewStart && t < viewEnd, m = 'pan';
        end
    end

    function tpp = ovTimePerPixel()
        p = getpixelposition(axOverview,true); tpp = max(dur,eps)/max(p(3),1);
    end

    function endOverviewDrag()
        set(fig,'WindowButtonMotionFcn',@(~,~)onHover(),'WindowButtonUpFcn','');
        ovMode = '';
        setStatus(sprintf('Viewing %s - %s.',fmtTime(viewStart),fmtTime(viewEnd)),'info');
    end

    % ---------------- typed time frame (From / To boxes) --------------------
    function beginFieldEdit(box)
        if isempty(vid), return; end
        if isEditing && ~isempty(typingBox) && ishghandle(typingBox) && ~isequal(typingBox,box)
            commitBox(typingBox);
        end
        typingBox = box; isEditing = true;
        % start with an empty box so you just type the number (old value is
        % kept in UserData and restored if you leave it empty)
        set(box,'UserData',get(box,'String'),'String','','Enable','on'); uicontrol(box);
        setStatus(['Type a time - e.g. 14 = 14 s, 1:05 = 1 min 5 s, 2m30 also works - ' ...
            'then press Enter or click Go.'],'info');
    end

    function syncTimeFields(a,b)
        % Show a range in the From / To boxes (skipped while you type in them).
        if isempty(edtFrom) || ~ishghandle(edtFrom), return; end
        if ~isequal(typingBox,edtFrom), set(edtFrom,'String',fmtTime(a)); end
        if ~isequal(typingBox,edtTo),   set(edtTo,'String',fmtTime(b));   end
    end

    function applyTimeFields(changed)
        % Reads the From / To boxes. Plain numbers are seconds (14 -> 00:14.000).
        % If the edited box makes the range impossible, the other end is moved
        % to keep the current window width instead of raising an error.
        if nargin < 1, changed = ''; end
        if isempty(vid), return; end
        a = parseTime(get(edtFrom,'String')); b = parseTime(get(edtTo,'String'));
        if isnan(a) || isnan(b)
            errordlg(sprintf(['That time could not be read.\n\nType seconds (e.g. 14 or 75.5), ' ...
                'minutes:seconds (e.g. 1:05), or 2m30.']),'Invalid time');
            refreshView(); return;
        end
        a = min(max(a,0),dur); b = min(max(b,0),dur);
        if b <= a
            w = max(viewEnd - viewStart, MINVIEW);
            if strcmp(changed,'tto'), a = max(0,b - w); else, b = min(dur,a + w); end
        end
        if b <= a
            errordlg(sprintf('The start time is at or past the end of the recording (%s).', ...
                fmtTime(dur)),'Invalid time frame');
            refreshView(); return;
        end
        setView(a,b);
        selectWholeView();
    end

    % ---------------- zoom / scroll / time frames ---------------------------
    function setView(a,b)
        if isempty(vid), return; end
        w = b - a;
        if w < MINVIEW, c=(a+b)/2; a=c-MINVIEW/2; b=c+MINVIEW/2; w=MINVIEW; end
        if w > dur, a=0; b=dur; end
        if a < 0,   b=b-a; a=0; end
        if b > dur, a=a-(b-dur); b=dur; end
        viewStart = max(0,a); viewEnd = min(dur,b); refreshView();
    end

    function zoomAbout(center,factor)
        if isempty(vid), return; end
        if isnan(center), center = (viewStart+viewEnd)/2; end
        w = (viewEnd-viewStart)*factor; w = max(MINVIEW,min(dur,w));
        frac = (center-viewStart)/max(viewEnd-viewStart,eps); frac = max(0,min(1,frac));
        a = center - frac*w; setView(a,a+w);
    end

    function selectWholeView()
        % Select the whole visible span (with its black sliders) so you can
        % label it straight away or trim it with the sliders.
        selStart = viewStart; selEnd = viewEnd; currentEvent = [];
        cursorTime = selStart;
        redrawEvents(); updateSelectionGraphics(); updateCursorGraphics(); showFrameAt(cursorTime);
        setStatus(sprintf(['Showing and selected %s - %s. Press R, B, P or T to label it, drag ' ...
            'the black sliders to trim it, or drag inside to pick a smaller part.'], ...
            fmtTime(selStart),fmtTime(selEnd)),'ok');
    end

    function fullView()
        setView(0,dur);
    end

    function zoomToSelection()
        if isnan(selStart) || isnan(selEnd) || selEnd <= selStart
            msgbox(sprintf(['There is no selection to zoom to.\n\nClick and drag across the ' ...
                'spectrogram or waveform first (or click an existing event).']), ...
                'Fit selection','help','replace');
            return;
        end
        pad = 0.1*(selEnd-selStart); setView(selStart-pad,selEnd+pad);
    end

    function scrollBy(frac)
        sh = frac*(viewEnd-viewStart); setView(viewStart+sh,viewEnd+sh);
    end

    function c = cursorCenter()
        if isnan(cursorTime), c=(viewStart+viewEnd)/2; else, c=cursorTime; end
    end

    function t = pointerTime()
        try
            cp = get(axSpec,'CurrentPoint'); t = min(max(cp(1,1),viewStart),viewEnd);
        catch
            t = cursorCenter();
        end
    end

    function goToTimeDialog()
        if isempty(vid), return; end
        answ = inputdlg({sprintf('Start time  (seconds, or m:ss)\nRecording length: %s', ...
            fmtTime(dur)),'End time  (seconds, or m:ss)'},'Go to time frame',[1 45], ...
            {fmtTime(viewStart),fmtTime(viewEnd)});
        if isempty(answ), return; end
        a = parseTime(answ{1}); b = parseTime(answ{2});
        if isnan(a) || isnan(b)
            errordlg(sprintf(['Those times could not be read.\n\nUse seconds (e.g. 75.5) or ' ...
                'minutes:seconds (e.g. 1:15.5).']),'Invalid time'); return;
        end
        if b <= a
            errordlg('The end time must be after the start time.','Invalid time frame'); return;
        end
        if a >= dur
            errordlg(sprintf('The start time is past the end of the recording (%s).', ...
                fmtTime(dur)),'Invalid time frame'); return;
        end
        b = min(b,dur);
        setView(a,b);
        selectWholeView();
    end

    function saveFrame()
        if isempty(vid), return; end
        def = sprintf('Frame %d',numel(frames)+1);
        answ = inputdlg({sprintf('Name for this time frame (%s - %s):', ...
            fmtTime(viewStart),fmtTime(viewEnd))},'Save time frame',[1 50],{def});
        if isempty(answ), return; end
        nm = strtrim(answ{1}); if isempty(nm), nm = def; end
        frames(end+1) = struct('name',nm,'a',viewStart,'b',viewEnd);
        updateFramesPopup(numel(frames)); dirty = true; updateTitle();
        setStatus(sprintf('Saved time frame "%s". Pick it from the drop-down to come back to it.',nm),'ok');
    end

    function onFramePick(src)
        v = get(src,'Value');
        if v <= 1 || v-1 > numel(frames), return; end
        f = frames(v-1);
        setView(f.a,f.b); cursorTime = f.a; updateCursorGraphics(); showFrameAt(f.a);
        setStatus(sprintf('Jumped to time frame "%s" (%s - %s).',f.name,fmtTime(f.a),fmtTime(f.b)),'info');
    end

    function removeFrame()
        v = get(popFrames,'Value');
        if v <= 1 || v-1 > numel(frames)
            msgbox('Pick a saved time frame from the drop-down first, then click "Remove frame".', ...
                'Remove time frame','help','replace'); return;
        end
        c = questdlg(sprintf('Remove the saved time frame "%s"?',frames(v-1).name), ...
            'Remove time frame','Remove','Cancel','Cancel');
        if ~strcmp(c,'Remove'), return; end
        frames(v-1) = []; updateFramesPopup(0); dirty = true; updateTitle();
        setStatus('Time frame removed.','info');
    end

    function updateFramesPopup(selFrame)
        items = {'Saved time frames...'};
        for q = 1:numel(frames)
            items{end+1} = sprintf('%s  [%s - %s]',frames(q).name, ...
                fmtTime(frames(q).a),fmtTime(frames(q).b)); %#ok<AGROW>
        end
        set(popFrames,'String',items,'Value',min(max(selFrame+1,1),numel(items)));
    end

    % ---------------- mouse on spectrogram / waveform / strips ---------------
    % Press on a black edge handle -> resize.  Press inside an event block ->
    % drag to move the whole block (a plain click selects it / types its text
    % on the strips).  Press elsewhere (or Shift + press) -> drag a new selection.
    function onAxDown(ax,panel)
        if nargin < 2, panel = ''; end
        if isempty(vid), return; end
        commitEdit();
        cp = get(ax,'CurrentPoint'); t = clampT(cp(1,1));
        shift = any(strcmp(get(fig,'CurrentModifier'),'shift'));
        dragPanel = panel; resizeUndoPushed = false; movePushed = false; moveIdx = [];
        dragMode = '';
        if ~shift, dragMode = hitEdge(ax,t); end
        if ~isempty(dragMode)
            if ~isempty(currentEvent), pushUndo(); resizeUndoPushed = true; end
        else
            hit = eventAtTime(t,~isempty(panel));
            if ~shift && ~isempty(hit)
                dragMode = 'move'; moveIdx = hit;
                moveOrig = [events(hit).start events(hit).end];
            else
                dragMode = 'select';
            end
        end
        dragAxes = ax; dragStartT = t; downPix = get(fig,'CurrentPoint'); didDrag = false;
        set(fig,'WindowButtonMotionFcn',@(~,~)cb(@onDrag),'WindowButtonUpFcn',@(~,~)cb(@onUp));
    end

    function onDrag()
        cp = get(dragAxes,'CurrentPoint'); t = clampT(cp(1,1));
        if norm(get(fig,'CurrentPoint')-downPix) > 3, didDrag = true; end
        if ~didDrag, return; end
        switch dragMode
            case 'left'
                selStart = max(0, min(t, selEnd-MINSEL));
                if ~isempty(currentEvent), events(currentEvent).start = selStart; redrawEvents(); end
            case 'right'
                selEnd = min(dur, max(t, selStart+MINSEL));
                if ~isempty(currentEvent), events(currentEvent).end = selEnd; redrawEvents(); end
            case 'move'
                if ~isempty(moveIdx) && ~movePushed
                    pushUndo(); movePushed = true;
                    currentEvent = moveIdx;
                end
                d = t - dragStartT; a = moveOrig(1)+d; b = moveOrig(2)+d;
                if a < 0,   b = b - a; a = 0; end
                if b > dur, a = a - (b-dur); b = dur; end
                selStart = a; selEnd = b;
                if ~isempty(moveIdx)
                    events(moveIdx).start = a; events(moveIdx).end = b; redrawEvents();
                end
            otherwise
                selStart = min(dragStartT,t); selEnd = max(dragStartT,t); currentEvent = [];
        end
        updateSelectionGraphics();
        if ~isnan(selStart), syncTimeFields(selStart,selEnd); end   % boxes follow the sliders
        drawnow limitrate;
    end

    function onUp()
        set(fig,'WindowButtonMotionFcn',@(~,~)onHover(),'WindowButtonUpFcn','');
        m0 = dragMode; dragMode = ''; panel = dragPanel; t = dragStartT;
        if any(strcmp(m0,{'left','right'}))
            if didDrag
                cursorTime = selStart;
                if ~isempty(currentEvent)
                    dirty = true; updateTitle(); updateEventList(); redrawEvents();
                    setStatus(sprintf('Resized %s.',eventLabel(currentEvent)),'ok');
                else
                    setStatus(sprintf('Selection %s - %s (%.3f s). Press R, B, P or T to label it.', ...
                        fmtTime(selStart),fmtTime(selEnd),selEnd-selStart),'info');
                end
            elseif resizeUndoPushed && ~isempty(undoStack)
                undoStack(end) = [];
            end
        elseif strcmp(m0,'move') && didDrag
            cursorTime = selStart;
            if ~isempty(moveIdx)
                dirty = true; updateTitle(); updateEventList(); redrawEvents();
                setStatus(sprintf('Moved %s (Ctrl+Z to undo).',eventLabel(moveIdx)),'ok');
            end
        elseif didDrag
            cursorTime = selStart; redrawEvents();
            setStatus(sprintf(['Selected %s - %s (%.3f s). Press R, B, P or T to label it, ' ...
                'Space to listen, or drag the black edges to adjust.'], ...
                fmtTime(selStart),fmtTime(selEnd),selEnd-selStart),'info');
        elseif ~isempty(panel)
            % plain click on the Transcript / Notes strip
            idx = moveIdx; if isempty(idx), idx = eventAtTime(t,true); end
            hasSel = ~isnan(selStart) && ~isnan(selEnd) && selEnd > selStart;
            if isempty(idx) && hasSel && isempty(currentEvent) && t >= selStart && t <= selEnd
                addEvent(TXT_TYPE,panel);          % text-only segment for the selection
            elseif ~isempty(idx)
                currentEvent = idx; selStart = events(idx).start; selEnd = events(idx).end;
                cursorTime = selStart; syncListSelection(); updateSelectionGraphics();
                startEditing(idx,panel);
            else
                setStatus(['To add text without an event: drag to select the words, ' ...
                    'then click this strip (or press T).'],'info');
            end
        else
            % plain click on the spectrogram / waveform
            idx = moveIdx; if isempty(idx), idx = eventAtTime(t,false); end
            if ~isempty(idx)
                currentEvent = idx; selStart = events(idx).start; selEnd = events(idx).end;
                cursorTime = selStart; syncListSelection();
                setStatus(sprintf(['Selected %s. Drag it to move it, drag its black edges to ' ...
                    'resize, R/B/P to change type, Delete to remove.'],eventLabel(idx)),'info');
            else
                selStart = NaN; selEnd = NaN; currentEvent = []; cursorTime = t;
            end
            redrawEvents();
        end
        moveIdx = [];
        updateSelectionGraphics(); updateCursorGraphics(); showFrameAt(cursorTime);
    end

    function edge = hitEdge(ax,t)
        edge = '';
        if isnan(selStart) || isnan(selEnd) || selEnd <= selStart, return; end
        tol = EDGE_PX*timePerPixel(ax);
        dL = abs(t-selStart); dR = abs(t-selEnd);
        if min(dL,dR) > tol, return; end
        if dL <= dR, edge = 'left'; else, edge = 'right'; end
    end

    function tpp = timePerPixel(ax)
        p = getpixelposition(ax,true); tpp = (viewEnd-viewStart)/max(p(3),1);
    end

    function onHover()
        % Show a resize pointer when hovering over a selection edge handle.
        try
            if isempty(vid) || isEditing || isempty(fig) || ~ishghandle(fig), return; end
            ptr = 'arrow';
            axs = [axSpec axWave axTrans axNotes];
            for q = 1:numel(axs)
                hx = axs(q);
                cp = get(hx,'CurrentPoint'); xl = get(hx,'XLim'); yl = get(hx,'YLim');
                if cp(1,1)>=xl(1) && cp(1,1)<=xl(2) && cp(1,2)>=yl(1) && cp(1,2)<=yl(2)
                    x = cp(1,1);
                    if ~isempty(hitEdge(hx,x))
                        ptr = 'left';
                    elseif ~isempty(eventAtTime(x,q > 2))
                        ptr = 'fleur';
                    end
                    break;
                end
            end
            cp = get(axOverview,'CurrentPoint');
            if cp(1,2) >= 0 && cp(1,2) <= 1 && cp(1,1) >= 0 && cp(1,1) <= dur
                m = ovHit(cp(1,1));
                if any(strcmp(m,{'left','right'})), ptr = 'left';
                elseif strcmp(m,'pan'), ptr = 'fleur'; end
            end
            if ~strcmp(get(fig,'Pointer'),ptr), set(fig,'Pointer',ptr); end
        catch
        end
    end

    function onResize()
        try
            if isempty(vid), return; end
            updateSelectionGraphics(); positionBoxes();
        catch
        end
    end

    % ---------------- keyboard ----------------------------------------------
    function onKey(e)
        k = e.Key;
        if strcmp(k,'f1'), showHelp(); return; end
        if isempty(vid)
            if any(strcmp(k,{'space','r','b','p','t'}))
                setStatus('Load a video first (click "Load video").','warn');
            end
            return;
        end
        if isEditing, return; end
        ctrl  = any(strcmp(e.Modifier,'control')) || any(strcmp(e.Modifier,'command'));
        shift = any(strcmp(e.Modifier,'shift'));
        if ctrl
            switch k
                case 'o', zoomAbout(cursorCenter(),2);
                case 'i', zoomAbout(cursorCenter(),0.5);
                case 'n', zoomToSelection();
                case 'a', fullView();
                case 'g', goToTimeDialog();
                case 's', saveAnnotations();
                case 'z', undo();
            end
            return;
        end
        switch k
            case 'space',                togglePlay();
            case 'escape'
                if isPlaying, stopPlay(); setStatus('Stopped.','info'); else, clearSelection(); end
            case {'delete','backspace'}, deleteCurrentEvent();
            case 'leftarrow',            stepCursor(-1,shift);
            case 'rightarrow',           stepCursor(1,shift);
            otherwise
                for q = 1:numel(eventTypes)
                    if ~isempty(eventTypes(q).key) && strcmp(k,eventTypes(q).key)
                        applyType(q); break;
                    end
                end
        end
    end

    function onScroll(e)
        if isempty(vid) || isEditing, return; end
        ctrl = any(ismember({'control','command'},get(fig,'CurrentModifier')));
        c = e.VerticalScrollCount;
        if ctrl, zoomAbout(pointerTime(),1.2^c); else, scrollBy(0.15*c); end
    end

    function stepCursor(direction,big)
        if isPlaying, return; end
        if big, step = 1; else, step = 1/max(vid.FrameRate,1); end
        t0 = cursorTime; if isnan(t0), t0 = viewStart; end
        cursorTime = clampT(t0 + direction*step);
        if cursorTime < viewStart || cursorTime > viewEnd
            w = viewEnd-viewStart; setView(cursorTime-w/2,cursorTime+w/2);
        end
        updateCursorGraphics(); showFrameAt(cursorTime);
        setStatus(sprintf('Cursor at %s',fmtTime(cursorTime)),'info');
    end

    function clearSelection()
        selStart = NaN; selEnd = NaN; currentEvent = [];
        redrawEvents(); updateSelectionGraphics();
        setStatus('Selection cleared.','info');
    end

    % ---------------- events ------------------------------------------------
    function applyType(q)
        if ~isempty(currentEvent) && (eventTypes(q).plain || ...
                eventTypes(events(currentEvent).type).plain)
            % adding text to an event's range, or an event inside a text
            % segment: keep the existing one and add a new one alongside
            currentEvent = [];
        end
        if ~isempty(currentEvent)
            k = currentEvent;
            if events(k).type == q
                setStatus(sprintf('%s is already a %s. Drag a new selection to add another event.', ...
                    eventLabel(k),eventTypes(q).name),'info');
                return;
            end
            pushUndo(); old = eventTypes(events(k).type).name;
            events(k).type = q; eventsChanged();
            setStatus(sprintf('Changed %s from %s to %s (Ctrl+Z to undo).', ...
                eventLabel(k),old,eventTypes(q).name),'ok');
        else
            addEvent(q);
        end
    end

    function addEvent(q,panel)
        if nargin < 2, panel = 'transcript'; end
        if isnan(selStart) || isnan(selEnd) || selEnd <= selStart
            msgbox(sprintf(['No time window is selected.\n\nClick and drag across the ' ...
                'spectrogram or waveform to select the stretch of speech you want to ' ...
                'label, then press %s again.'],upper(eventTypes(q).key)), ...
                'Select a section first','help','replace');
            return;
        end
        if selEnd - selStart < 0.02
            c = questdlg(sprintf(['This selection is only %.0f ms long - it may have been ' ...
                'an accidental click-drag.\n\nAdd a %s event anyway?'], ...
                1000*(selEnd-selStart),eventTypes(q).name),'Very short event','Add','Cancel','Add');
            if ~strcmp(c,'Add'), return; end
        end
        ov = overlappingEvents(selStart,selEnd);
        if ~isempty(ov)   % text-only segments never count as clashes
            ov = ov(~[eventTypes([events(ov).type]).plain]);
        end
        if eventTypes(q).plain, ov = []; end
        if ~isempty(ov)
            c = questdlg(sprintf(['This selection overlaps %d existing event(s), e.g. %s.\n\n' ...
                'Overlapping labels are allowed (for example a repetition inside a block).\n' ...
                'Add the new %s event anyway?'],numel(ov),eventLabel(ov(1)),eventTypes(q).name), ...
                'Overlapping events','Add anyway','Cancel','Add anyway');
            if ~strcmp(c,'Add anyway'), return; end
        end
        pushUndo();
        s = struct('start',selStart,'end',selEnd,'type',q,'transcript','','notes','');
        events(end+1) = s; currentEvent = numel(events);
        eventsChanged(); updateSelectionGraphics();
        startEditing(currentEvent,panel);
    end

    function eventsChanged()
        dirty = true; updateTitle(); updateEventList(); redrawEvents();
    end

    function idx = eventAtTime(t,includePlain)
        % Shortest event containing t (so nested events stay clickable).
        idx = [];
        if isempty(events), return; end
        if nargin < 2, includePlain = true; end
        hit = find([events.start] <= t & [events.end] >= t);
        if ~includePlain && ~isempty(hit)
            hit = hit(~[eventTypes([events(hit).type]).plain]);
        end
        if isempty(hit), return; end
        [~,m] = min([events(hit).end] - [events(hit).start]);
        idx = hit(m);
    end

    function idx = overlappingEvents(a,b)
        idx = [];
        if isempty(events), return; end
        idx = find([events.start] < b & [events.end] > a);
    end

    function deleteCurrentEvent()
        if isempty(currentEvent)
            msgbox(sprintf(['No event is selected.\n\nClick an event block on the timeline ' ...
                '(or pick it in the Events list), then press Delete.']), ...
                'Delete event','help','replace');
            return;
        end
        k = currentEvent; tr = events(k).transcript; if isempty(tr), tr = '(none)'; end
        c = questdlg(sprintf('Delete %s?\n\nTranscript: %s\n\nYou can undo this with Ctrl+Z.', ...
            eventLabel(k),tr),'Delete event','Delete','Cancel','Cancel');
        if ~strcmp(c,'Delete'), return; end
        commitEditQuiet(); pushUndo();
        events(k) = []; currentEvent = []; selStart = NaN; selEnd = NaN;
        eventsChanged(); updateSelectionGraphics();
        setStatus('Event deleted (Ctrl+Z to undo).','ok');
    end

    function changeTypeDialog()
        if isempty(currentEvent)
            msgbox('Select an event first (click it on the timeline or in the Events list).', ...
                'Change type','help','replace'); return;
        end
        k = currentEvent;
        [sel,ok] = listdlg('PromptString','Choose the new event type:', ...
            'SelectionMode','single','ListString',{eventTypes.name}, ...
            'InitialValue',events(k).type,'Name','Change event type','ListSize',[240 120]);
        if ok, applyType(sel); end
    end

    function editCurrentEvent()
        if isempty(currentEvent)
            msgbox('Select an event first (click it on the timeline or in the Events list).', ...
                'Edit text','help','replace'); return;
        end
        k = currentEvent;
        if events(k).start < viewStart || events(k).end > viewEnd, selectEvent(k,true); end
        startEditing(k,'transcript');
    end

    function playCurrentEvent()
        if isempty(currentEvent)
            msgbox('Select an event first (click it on the timeline or in the Events list).', ...
                'Play event','help','replace'); return;
        end
        stopPlay(); selectEvent(currentEvent,true); startPlay();
    end

    function selectEvent(k,ensureVisible)
        currentEvent = k; selStart = events(k).start; selEnd = events(k).end;
        cursorTime = selStart;
        if ensureVisible && (selStart < viewStart || selEnd > viewEnd)
            w = max(viewEnd-viewStart,1.5*(selEnd-selStart)); c = (selStart+selEnd)/2;
            setView(c-w/2,c+w/2);
        else
            redrawEvents(); updateSelectionGraphics(); updateCursorGraphics();
        end
        showFrameAt(cursorTime); syncListSelection();
    end

    % ---------------- events list -------------------------------------------
    function updateEventList()
        if isempty(events)
            listMap = [];
            set(lstEvents,'Data',cell(0,5),'BackgroundColor',[1 1 1]); return;
        end
        [~,ord] = sort([events.start]); listMap = ord;
        n = numel(ord); data = cell(n,5); bgc = ones(n,3);
        for q = 1:n
            ev = events(ord(q)); et = eventTypes(ev.type);
            tr = strrep(ev.transcript,sprintf('\n'),' ');
            num = sprintf('%d',q);
            if ~isempty(currentEvent) && ord(q) == currentEvent, num = ['> ' num]; end
            if et.plain, tyName = 'text'; else, tyName = shortType(ev.type); end
            data(q,:) = {num, fmtTime(ev.start), sprintf('%.2f',ev.end-ev.start), tyName, tr};
            if ~et.plain, bgc(q,:) = 0.65 + 0.35*et.color; end   % light tint
        end
        set(lstEvents,'Data',data,'BackgroundColor',bgc);
    end

    function syncListSelection()
        if isempty(listMap), return; end
        updateEventList();          % moves the '>' marker to the selected row
    end

    function onListSelect(~,e)
        if isempty(vid) || isempty(listMap) || isempty(e.Indices), return; end
        commitEdit();
        r = e.Indices(1,1); if r < 1 || r > numel(listMap), return; end
        k = listMap(r);
        selectEvent(k,true);
        setStatus(sprintf('Selected %s. Press Space or "Play event" to hear it.',eventLabel(k)),'info');
    end

    function s = eventLabel(k)
        pos = find(listMap == k,1); if isempty(pos), pos = k; end
        s = sprintf('event #%d (%s, %s - %s)',pos,eventTypes(events(k).type).name, ...
            fmtTime(events(k).start),fmtTime(events(k).end));
    end

    function s = eventLabelShort(k)
        pos = find(listMap == k,1); if isempty(pos), pos = k; end
        s = sprintf('#%d %s',pos,eventTypes(events(k).type).name);
    end

    function s = shortType(ti)
        nm = eventTypes(ti).name; s = upper(nm(1:min(4,numel(nm))));
    end

    % ---------------- annotation text editing -------------------------------
    % When an event is selected, white text boxes sit on top of its block in
    % the Transcript and Notes strips. Click one and type; press Enter, Tab,
    % or click anywhere else to save.
    function startEditing(idx,panel)
        if isempty(vid) || isempty(idx) || idx > numel(events), return; end
        if ~isequal(currentEvent,idx)
            currentEvent = idx; selStart = events(idx).start; selEnd = events(idx).end;
            syncListSelection(); updateSelectionGraphics();
        end
        redrawEvents();                 % positions the boxes
        beginTyping(panel);
    end

    function beginTyping(panel)
        if isempty(vid) || isempty(currentEvent), return; end
        if isEditing && ~isempty(typingBox) && ishghandle(typingBox)
            commitBox(typingBox);
        end
        if strcmp(panel,'notes'), box = notesBox; fld = 'notes';
        else,                     box = transBox; fld = 'transcript'; end
        k = currentEvent;
        set(box,'String',events(k).(fld),'UserData',k);
        typingBox = box; isEditing = true;
        positionBoxes();
        if ~strcmp(get(box,'Visible'),'on'), typingBox = []; isEditing = false; return; end
        set(box,'Enable','on');
        uicontrol(box);                 % put the text cursor in the box
        setStatus(sprintf(['Typing the %s for %s - press Enter (or click anywhere ' ...
            'else) to save.'],panel,eventLabel(k)),'info');
    end

    function commitBox(box)
        if isempty(box) || ~ishghandle(box), return; end
        if any(strcmp(get(box,'Tag'),{'tfrom','tto'}))
            set(box,'Enable','inactive');
            if isequal(typingBox,box), typingBox = []; isEditing = false; end
            if isempty(strtrim(char(get(box,'String'))))   % left empty: restore
                set(box,'String',get(box,'UserData'));
            end
            applyTimeFields(get(box,'Tag')); return;
        end
        idx = get(box,'UserData'); field = get(box,'Tag');
        str = get(box,'String');
        if iscell(str), str = strjoin(str,' ');
        elseif size(str,1) > 1, str = strjoin(cellstr(str),' '); end
        str = strtrim(str);
        set(box,'Enable','inactive','Visible','off');
        if isequal(typingBox,box), typingBox = []; isEditing = false; end
        if isempty(idx) || idx > numel(events), return; end
        if ~strcmp(events(idx).(field),str)
            pushUndo();
            events(idx).(field) = str;
            dirty = true; updateTitle(); updateEventList(); redrawEvents();
            setStatus(sprintf('Saved %s for %s.',field,eventLabel(idx)),'ok');
        end
    end

    function commitEdit()
        % Let a pending edit-box callback run first (macOS commits the text
        % on focus loss), then save whatever is still being typed.
        drawnow;
        if isEditing && ~isempty(typingBox) && ishghandle(typingBox)
            commitBox(typingBox);
        end
        isEditing = false; typingBox = [];
    end

    function commitEditQuiet()
        isEditing = false; typingBox = [];
        if ~isempty(transBox) && ishghandle(transBox), set(transBox,'Enable','inactive'); end
        if ~isempty(notesBox) && ishghandle(notesBox), set(notesBox,'Enable','inactive'); end
        if ~isempty(edtFrom) && ishghandle(edtFrom), set(edtFrom,'Enable','inactive'); end
        if ~isempty(edtTo) && ishghandle(edtTo), set(edtTo,'Enable','inactive'); end
    end

    function positionBoxes()
        % A white text box is shown only while you type in it.
        boxes = {transBox, axTrans; notesBox, axNotes};
        for q = 1:2
            box = boxes{q,1};
            if isempty(box) || ~ishghandle(box), continue; end
            if ~isequal(typingBox,box), set(box,'Visible','off'); continue; end
            k = get(box,'UserData');
            show = ~isempty(vid) && ~isempty(k) && k <= numel(events) && ...
                events(k).end > viewStart && events(k).start < viewEnd;
            if ~show, commitBox(box); set(box,'Visible','off'); continue; end
            x0 = max(events(k).start,viewStart); x1 = min(events(k).end,viewEnd);
            set(box,'Position',dataRangeToPix(boxes{q,2},x0,x1),'Visible','on');
        end
    end

    function pos = dataRangeToPix(ax,x0,x1)
        axpix = getpixelposition(ax,true); xl = get(ax,'XLim');
        f0 = (x0-xl(1))/diff(xl); f1 = (x1-xl(1))/diff(xl);
        f0 = max(0,min(1,f0)); f1 = max(0,min(1,f1));
        pw = max(200,(f1-f0)*axpix(3)); pw = min(pw,axpix(3));
        px = axpix(1)+f0*axpix(3);
        px = min(px, axpix(1)+axpix(3)-pw);          % keep inside the strip
        ph = min(28, axpix(4)-6);
        pos = [px axpix(2)+(axpix(4)-ph)/2 pw ph];
    end

    % ---------------- undo --------------------------------------------------
    function pushUndo()
        undoStack{end+1} = events;
        if numel(undoStack) > MAXUNDO, undoStack(1) = []; end
    end

    function undo()
        if isempty(undoStack), setStatus('Nothing to undo.','warn'); return; end
        commitEditQuiet();
        events = undoStack{end}; undoStack(end) = [];
        currentEvent = []; selStart = NaN; selEnd = NaN;
        eventsChanged(); updateSelectionGraphics();
        setStatus('Undid the last change.','ok');
    end

    % ---------------- playback ----------------------------------------------
    % One audioplayer is built per video and reused (play(player,[s0 s1])).
    % Creating a new player on every key press could leave macOS playing
    % silently; reusing one device stream is more reliable.
    function buildPlayer()
        try, if ~isempty(mainPlayer), stop(mainPlayer); end; catch, end
        mainPlayer = [];
        try
            if audiodevinfo(0) < 1
                errordlg(['MATLAB cannot find an audio output device, so playback ' ...
                    'will be silent. Check your Mac''s sound output and restart MATLAB.'], ...
                    'No audio output');
            end
        catch
        end
        try
            mainPlayer = audioplayer(audio,fs);
        catch err
            errordlg(sprintf('Could not set up audio playback:\n\n%s',err.message),'Audio error');
        end
    end

    function togglePlay()
        if isPlaying, pausePlay(); else, startPlay(); end
    end

    function startPlay()
        if isempty(vid) || isPlaying, return; end
        if isempty(mainPlayer), buildPlayer(); end
        if isempty(mainPlayer), return; end
        hasSel = ~isnan(selStart) && ~isnan(selEnd) && selEnd > selStart;
        if hasSel && ~isnan(cursorTime) && cursorTime > selStart && cursorTime < selEnd - 0.01
            t0 = cursorTime; t1 = selEnd; whatTxt = 'rest of selection';   % resume
        elseif hasSel
            t0 = selStart; t1 = selEnd; whatTxt = 'selection';
        elseif ~isnan(cursorTime)
            t0 = cursorTime; t1 = dur; whatTxt = 'from cursor';
        else
            t0 = viewStart; t1 = dur; whatTxt = 'from view start';
        end
        s0 = max(1,round(t0*fs)+1); s1 = min(numel(audio),round(t1*fs));
        if s1 <= s0, return; end
        try
            set(mainPlayer,'SampleRate',round(fs*playSpeed));
        catch
            playSpeed = 1; set(popSpeed,'Value',1);
            try, set(mainPlayer,'SampleRate',fs); catch, end
            setStatus('This audio device does not support that speed - playing at 1x.','warn');
        end
        playStartTime = 0; playEndTime = t1; isPlaying = true;
        deleteValid(playGfx);
        pl1 = plot(axSpec,[t0 t0],get(axSpec,'YLim'),'-','Color',PLAY_COLOR,'LineWidth',1.5);
        pl2 = plot(axWave,[t0 t0],[-1 1],'-','Color',PLAY_COLOR,'LineWidth',1.5);
        playGfx = [pl1 pl2]; set(playGfx,'PickableParts','none','HitTest','off');
        try
            vid.CurrentTime = max(0, min(t0, vid.Duration - 1/max(vid.FrameRate,1)));
        catch
        end
        try
            play(mainPlayer,[s0 s1]);
        catch err
            isPlaying = false; deleteValid(playGfx);
            errordlg(sprintf('Audio playback failed:\n\n%s',err.message),'Audio error'); return;
        end
        start(playTimer);
        msg = sprintf('Playing %s at %.2fx...   Space = pause.',whatTxt,playSpeed);
        if max(abs(audio(s0:s1))) < 1e-4, msg = [msg '  (This part of the recording is silent.)']; end
        setStatus(msg,'info');
    end

    function playTick()
        if ~isPlaying, return; end
        player = mainPlayer;
        if isempty(player) || ~isplaying(player), stopPlay(); return; end
        t = currentPlayTime(player);                     % audio = master clock
        if t >= playEndTime, stopPlay(); return; end
        if t > viewEnd || t < viewStart
            w = viewEnd-viewStart; viewStart = max(0,t-0.1*w); viewEnd = viewStart+w; refreshView();
        end
        if ~isempty(playGfx) && all(ishghandle(playGfx))
            set(playGfx(1),'XData',[t t]); set(playGfx(2),'XData',[t t]);
        end
        if ~isempty(ovCur) && all(isgraphics(ovCur)), set(ovCur,'XData',[t t]); end
        advanceFrameTo(t); drawnow limitrate;
    end

    function t = currentPlayTime(player)
        % CurrentSample is the absolute sample index in the full recording.
        t = playStartTime + double(player.CurrentSample - 1)/fs;
    end

    function pausePlay()
        % Pause and put the cursor exactly where playback stopped; the next
        % Space press resumes from there.
        if ~isPlaying, return; end
        t = NaN;
        try, t = currentPlayTime(mainPlayer); catch, end
        stopPlay();
        if isfinite(t)
            cursorTime = clampT(t);
            updateCursorGraphics(); showFrameAt(cursorTime);
        end
        setStatus(sprintf('Paused at %s. Press Space to continue from here.',fmtTime(cursorTime)),'info');
    end

    function stopPlay()
        try, stop(playTimer); end %#ok<TRYNC>
        try, if ~isempty(mainPlayer), stop(mainPlayer); end; catch, end
        isPlaying = false; deleteValid(playGfx); playGfx = gobjects(0);
        if ~isempty(ovCur) && all(isgraphics(ovCur)), set(ovCur,'XData',[cursorTime cursorTime]); end
    end

    function onSpeed(src)
        playSpeed = SPEEDS(get(src,'Value'));
        msg = sprintf('Playback speed set to %.2fx (applies the next time you press Space).',playSpeed);
        if playSpeed < 1, msg = [msg ' Slower speeds also lower the pitch.']; end
        setStatus(msg,'info');
    end

    function advanceFrameTo(t)
        % Playback frame update: decode forward instead of seeking every tick.
        if isempty(vid) || isempty(imgVideo) || ~ishghandle(imgVideo), return; end
        fp = 1/max(vid.FrameRate,1);
        try
            if t < vid.CurrentTime - fp || t > vid.CurrentTime + 0.5
                vid.CurrentTime = max(0, min(t, vid.Duration - fp));
            end
            fr = [];
            while hasFrame(vid) && vid.CurrentTime <= t
                fr = readFrame(vid);
            end
            if ~isempty(fr), set(imgVideo,'CData',fr); end
        catch
        end
    end

    function showFrameAt(t)
        % Random-access frame display (clicks, drags, pauses).
        if isempty(vid) || isempty(imgVideo) || ~ishghandle(imgVideo) || isnan(t), return; end
        tt = max(0,min(t, vid.Duration - 1/max(vid.FrameRate,1)));
        try
            vid.CurrentTime = tt; set(imgVideo,'CData',readFrame(vid));
        catch
        end
    end

    % ---------------- save / load annotations -------------------------------
    function ok = saveAnnotations()
        ok = false;
        if isempty(vid), errordlg('Load a video first.','Save'); return; end
        commitEdit();
        if isempty(events)
            c = questdlg('There are no events to save yet. Save an empty annotation file anyway?', ...
                'Nothing to save','Save anyway','Cancel','Cancel');
            if ~strcmp(c,'Save anyway'), return; end
        end
        [vdir,base] = fileparts(videoPath);
        defName = fullfile(vdir,[base '_annot-disfluencies.xlsx']);
        [fn,fp] = uiputfile({'*.xlsx','Excel workbook (*.xlsx)'},'Save annotations',defName);
        if isequal(fn,0), return; end
        out = fullfile(fp,fn);
        if isempty(events)
            starts = zeros(0,1); ends = zeros(0,1); durs = zeros(0,1);
            types = cell(0,1); trans = cell(0,1); notesv = cell(0,1);
        else
            [~,ord] = sort([events.start]); ev = events(ord);
            starts = [ev.start]'; ends = [ev.end]'; durs = ends - starts;
            types  = arrayfun(@(x)eventTypes(x.type).name,ev,'UniformOutput',false)';
            trans  = arrayfun(@(x)x.transcript,ev,'UniformOutput',false)';
            notesv = arrayfun(@(x)x.notes,ev,'UniformOutput',false)';
        end
        T = table(starts,ends,types,trans,notesv,durs, ...
            'VariableNames',{'starts','ends','event_type','transcript','notes','duration_s'});
        try
            if exist(out,'file'), delete(out); end      % avoid stale sheets
            writetable(T,out,'Sheet','events');
            if ~isempty(frames)
                F = table({frames.name}',[frames.a]',[frames.b]', ...
                    'VariableNames',{'name','start_s','end_s'});
                writetable(F,out,'Sheet','time_frames');
            end
        catch err
            errordlg(sprintf(['Could not save the annotations:\n\n%s\n\nIf the file is open ' ...
                'in Excel, close it and try again.'],err.message),'Save error');
            return;
        end
        dirty = false; updateTitle(); ok = true;
        setStatus(sprintf('Saved %d event(s) to %s',numel(events),fn),'ok');
        msgbox(sprintf('Saved %d event(s) and %d time frame(s) to:\n\n%s', ...
            numel(events),numel(frames),out),'Annotations saved','help','replace');
    end

    function loadAnnotations()
        if isempty(vid)
            errordlg('Load the matching video first, then load its annotations.','Load annotations');
            return;
        end
        commitEdit();
        [fn,fp] = uigetfile({'*.xlsx;*.xls;*.csv','Annotation tables (*.xlsx, *.xls, *.csv)'}, ...
            'Load annotations',fileparts(videoPath));
        if isequal(fn,0), return; end
        f = fullfile(fp,fn);
        try
            T = readtable(f);
        catch err
            errordlg(sprintf('Could not read the file:\n\n%s',err.message),'Load error'); return;
        end
        req = {'starts','ends','event_type','transcript','notes'};
        vn = lower(T.Properties.VariableNames); col = zeros(1,5);
        for j = 1:5, m = find(strcmp(vn,req{j}),1); if ~isempty(m), col(j) = m; end; end
        if any(col == 0)
            if width(T) >= 5
                col = 1:5;
            else
                errordlg(sprintf(['This file does not look like an annotation file.\n\n' ...
                    'Expected columns: starts, ends, event_type, transcript, notes.']), ...
                    'Invalid file'); return;
            end
        end

        typesBackup = eventTypes;
        newEv = struct('start',{},'end',{},'type',{},'transcript',{},'notes',{});
        nBad = 0; nOut = 0;
        for i = 1:height(T)
            st = toNum(T{i,col(1)}); en = toNum(T{i,col(2)});
            if isnan(st) || isnan(en) || en <= st, nBad = nBad + 1; continue; end
            if st >= dur, nOut = nOut + 1; continue; end
            en = min(en,dur);
            nm = toStr(T{i,col(3)}); if isempty(nm), nm = 'unlabelled'; end
            tr = toStr(T{i,col(4)}); nt = toStr(T{i,col(5)});
            ti = find(strcmpi({eventTypes.name},nm),1);
            if isempty(ti)
                eventTypes(end+1) = struct('name',nm,'key','','style','-.','color',[0.5 0.5 0.5],'plain',false); %#ok<AGROW>
                ti = numel(eventTypes);
            end
            newEv(end+1) = struct('start',st,'end',en,'type',ti,'transcript',tr,'notes',nt); %#ok<AGROW>
        end
        skipNote = '';
        if nBad > 0, skipNote = [skipNote sprintf('\n - %d row(s) with missing or invalid times will be skipped.',nBad)]; end
        if nOut > 0, skipNote = [skipNote sprintf('\n - %d row(s) start after the end of this video and will be skipped.',nOut)]; end
        if isempty(newEv)
            eventTypes = typesBackup;
            errordlg(sprintf('No usable events were found in "%s".%s',fn,skipNote),'Nothing to load');
            return;
        end

        mergeMode = false;
        if ~isempty(events)
            c = questdlg(sprintf(['This video already has %d annotation(s).\n\n"%s" contains ' ...
                '%d event(s). Do you want to add them to your current annotations (merge), ' ...
                'or replace your current annotations with them?'],numel(events),fn,numel(newEv)), ...
                'Load annotations','Merge with current','Replace current','Cancel','Merge with current');
            switch c
                case 'Merge with current'
                    mergeMode = true;
                    [nOverlap,isDup] = compareToCurrent(newEv);
                    msg = sprintf(['Are you sure you want to continue?\n\n' ...
                        'The annotations you are uploading may clash with the ones already on ' ...
                        'this video:\n' ...
                        ' - %d imported event(s) overlap in time with existing events\n' ...
                        ' - %d exact duplicate(s) will be skipped\n' ...
                        ' - %d new event(s) will be added%s\n\n' ...
                        'Overlapping events are kept side by side - nothing already on the ' ...
                        'timeline is changed or deleted. You can undo the whole import with Ctrl+Z.'], ...
                        nOverlap,sum(isDup),sum(~isDup),skipNote);
                    c2 = questdlg(msg,'Possible clash with current annotations', ...
                        'Continue','Cancel','Cancel');
                    if ~strcmp(c2,'Continue'), eventTypes = typesBackup; return; end
                    newEv = newEv(~isDup);
                    if isempty(newEv)
                        eventTypes = typesBackup;
                        msgbox('Every event in that file is already on the timeline - nothing was added.', ...
                            'Nothing new','help','replace');
                        return;
                    end
                case 'Replace current'
                    c2 = questdlg(sprintf(['Replace all %d current annotation(s) with the %d ' ...
                        'event(s) from "%s"?%s\n\nYou can undo this with Ctrl+Z.'], ...
                        numel(events),numel(newEv),fn,skipNote),'Replace annotations', ...
                        'Replace','Cancel','Cancel');
                    if ~strcmp(c2,'Replace'), eventTypes = typesBackup; return; end
                otherwise
                    eventTypes = typesBackup; return;
            end
        end

        commitEditQuiet(); pushUndo();
        if mergeMode
            events = [events newEv]; dirty = true; verb = 'Merged';
        else
            events = newEv; dirty = false; verb = 'Loaded';
        end

        % saved time frames (xlsx only; silently skipped if absent)
        nFr = 0;
        try
            Fr = readtable(f,'Sheet','time_frames');
            for i = 1:height(Fr)
                nm = toStr(Fr{i,1}); a = toNum(Fr{i,2}); b = toNum(Fr{i,3});
                if isempty(nm) || isnan(a) || isnan(b) || b <= a, continue; end
                if any(strcmp({frames.name},nm)), continue; end
                frames(end+1) = struct('name',nm,'a',a,'b',min(b,dur)); %#ok<AGROW>
                nFr = nFr + 1;
            end
        catch
        end

        currentEvent = []; selStart = NaN; selEnd = NaN;
        updateEventList(); updateFramesPopup(0); updateTitle(); refreshView();
        setStatus(sprintf('%s %d event(s) from %s.',verb,numel(newEv),fn),'ok');
        extra = '';
        if nFr > 0, extra = sprintf('\n%d saved time frame(s) were also added.',nFr); end
        if ~mergeMode && ~isempty(skipNote), extra = [extra sprintf('\n\nNote:%s',skipNote)]; end
        msgbox(sprintf('%s %d event(s) from "%s".%s',verb,numel(newEv),fn,extra), ...
            'Annotations loaded','help','replace');
    end

    function [nOverlap,isDup] = compareToCurrent(newEv)
        nOverlap = 0; isDup = false(1,numel(newEv));
        es = [events.start]; ee = [events.end]; et = [events.type];
        for i = 1:numel(newEv)
            d = abs(es-newEv(i).start) < 1e-3 & abs(ee-newEv(i).end) < 1e-3 & et == newEv(i).type;
            if any(d), isDup(i) = true; continue; end
            if any(es < newEv(i).end & ee > newEv(i).start), nOverlap = nOverlap + 1; end
        end
    end

    function ok = confirmUnsaved(actionText)
        ok = true;
        if ~dirty, return; end
        c = questdlg(sprintf('You have unsaved annotation changes.\n\nSave them before %s?', ...
            actionText),'Unsaved changes','Save','Don''t save','Cancel','Save');
        switch c
            case 'Save',       ok = saveAnnotations();
            case 'Don''t save', ok = true;
            otherwise,         ok = false;
        end
    end

    % ---------------- help windows ------------------------------------------
    function showHelp()
        msg = sprintf([ ...
            'GETTING STARTED\n' ...
            '  1. Click "Load video" and pick a recording (.mp4, .mov, .avi ...).\n' ...
            '  2. Click and drag across the spectrogram or waveform to select speech.\n' ...
            '  3. Press R (repetition), B (block) or P (prolongation) to label it.\n' ...
            '  4. Type what was said in the white box on the Transcript strip, and any\n' ...
            '     notes in the box on the Notes strip. Enter or clicking elsewhere saves.\n\n' ...
            'TRANSCRIPT / NOTES WITHOUT AN EVENT\n' ...
            '  - Select the words, then press T or click the Transcript / Notes strip.\n' ...
            '    This makes a text-only segment: grey outline, no coloured block, and\n' ...
            '    it does not count as a disfluency. Disfluency events can sit inside it.\n\n' ...
            'ADJUSTING AN EVENT\n' ...
            '  - Click an event block (or pick it in the Events list) to select it.\n' ...
            '  - Drag the black handles on its left / right edge to shorten or lengthen it\n' ...
            '    (on the spectrogram, waveform, Transcript or Notes strip).\n' ...
            '  - Drag the block itself to move it. Shift + drag makes a new selection\n' ...
            '    on top of an existing event instead.\n' ...
            '  - Press R / B / P while it is selected to change its type.\n' ...
            '  - Press Delete to remove it. Ctrl+Z undoes any change.\n\n' ...
            'LOOKING AT A SPECIFIC TIME FRAME\n' ...
            '  - "Go to time..." (Ctrl+G): type a start and end time (75.5 or 1:15.5).\n' ...
            '  - Zoom in / Zoom out, "Fit selection", "Full recording", or Ctrl + wheel.\n' ...
            '  - Click the Overview strip to jump; drag the blue window to slide it, or\n' ...
            '    drag its black edge handles to make the visible window wider / narrower.\n' ...
            '  - Type exact times in the "Show from ... to ..." boxes above the Overview\n' ...
            '    and press Enter or Go.\n' ...
            '  - "Save frame" stores the current view; pick it from the drop-down to\n' ...
            '    return to it. Saved frames are written into the annotation file.\n\n' ...
            'PLAYBACK\n' ...
            '  - Space plays the selection (or from the cursor); Space again pauses and\n' ...
            '    leaves the cursor exactly where it stopped, so Space resumes from there.\n' ...
            '  - Left / Right arrows step one video frame; Shift + arrows step 1 second.\n' ...
            '  - The speed menu gives slow-motion review (slower speeds lower the pitch).\n\n' ...
            'SAVING AND LOADING\n' ...
            '  - "Save annotations" (Ctrl+S) writes an Excel file next to the video.\n' ...
            '  - "Load annotations" can MERGE a file into the current annotations or\n' ...
            '    REPLACE them. You are warned about clashes before anything changes.\n' ...
            '  - You are reminded to save before closing or loading another video.\n\n' ...
            'READING THE DISPLAY\n' ...
            '  - Blocks are colour-coded: red = repetition, blue = block,\n' ...
            '    green = prolongation (outline: solid / dashed / dotted). All text is black.\n' ...
            '  - The grey shading with black handles is the current selection.\n' ...
            '  - Event times come from the audio sample clock, not video frames.']);
        msgbox(msg,[APPNAME ' - Help'],'help','replace');
    end

    function showShortcuts()
        % Custom two-column shortcuts window, sized to its content.
        old = findall(0,'Tag','DA_Shortcuts');
        if ~isempty(old), figure(old(1)); return; end

        rows = { ...
            'H','PLAYBACK',[]; ...
            'K','Space','Play the selection (or from the cursor)'; ...
            'K','Space again','Pause - cursor stays exactly where it stopped'; ...
            'K','Esc','Stop, or clear the selection when stopped'; ...
            'K','Left / Right','Step one video frame'; ...
            'K','Shift + Left / Right','Step one second'; ...
            'H','ZOOM & NAVIGATE',[]; ...
            'K','Ctrl/Cmd + I  /  O','Zoom in / out'; ...
            'K','Ctrl/Cmd + wheel','Zoom around the mouse pointer'; ...
            'K','Mouse wheel','Scroll through time'; ...
            'K','Ctrl/Cmd + N','Zoom to the selection'; ...
            'K','Ctrl/Cmd + A','Show the whole recording'; ...
            'K','Ctrl/Cmd + G','Go to a specific time frame'; ...
            'K','Click Overview strip','Jump to that point'; ...
            'K','Drag Overview edges','Widen / narrow the visible window'; ...
            'K','Drag Overview window','Slide the view through the recording'; ...
            'K','From / To boxes','Type a time frame, then Enter or Go'; ...
            'H','SELECTING & EDITING',[]; ...
            'K','Click + drag','Select a time window'; ...
            'K','Drag a black handle','Resize the selection / event'; ...
            'K','Drag an event block','Move the whole block'; ...
            'K','Shift + drag','New selection on top of an event'; ...
            'K','Click an event','Select it'; ...
            'K','Click Transcript / Notes','Type text for that event / selection'; ...
            'K','Delete','Delete the selected event'; ...
            'K','Ctrl/Cmd + Z','Undo'; ...
            'K','Ctrl/Cmd + S','Save annotations'; ...
            'K','F1','Full help'; ...
            'H','EVENT KEYS',[]};
        for q = 1:numel(eventTypes)
            if isempty(eventTypes(q).key), kd = '(no key)'; else, kd = upper(eventTypes(q).key); end
            rows(end+1,:) = {'E',kd,q}; %#ok<AGROW>
        end

        lineH = 22; gapH = 10; figW = 560; colKey = 40; colDesc = 230;
        totalH = 16;
        for q = 1:size(rows,1)
            if strcmp(rows{q,1},'H') && q > 1, totalH = totalH + gapH; end
            totalH = totalH + lineH;
        end
        btnH = 54; figHt = totalH + btnH + 10;

        mp = getpixelposition(fig);
        fx = mp(1) + (mp(3)-figW)/2; fy = mp(2) + (mp(4)-figHt)/2;
        sf = figure('Name','Keyboard & mouse shortcuts','NumberTitle','off', ...
            'MenuBar','none','ToolBar','none','Resize','off','Color',[1 1 1], ...
            'Units','pixels','Position',[fx fy figW figHt],'Tag','DA_Shortcuts', ...
            'HandleVisibility','off','KeyPressFcn',@closeOnKey);
        ax = axes('Parent',sf,'Units','pixels','Position',[0 btnH figW totalH+10], ...
            'XLim',[0 figW],'YLim',[0 totalH+10],'YDir','reverse','Visible','off', ...
            'HandleVisibility','off');

        y = 14;
        for q = 1:size(rows,1)
            switch rows{q,1}
                case 'H'
                    if q > 1, y = y + gapH; end
                    text(ax,24,y,rows{q,2},'FontWeight','bold','FontSize',11, ...
                        'Color',[0 0 0],'VerticalAlignment','top');
                    line(ax,[24 figW-24],[y+lineH-3 y+lineH-3],'Color',[0.8 0.8 0.8]);
                case 'K'
                    text(ax,colKey,y,rows{q,2},'FontWeight','bold','FontSize',10, ...
                        'Color',[0 0 0],'VerticalAlignment','top');
                    text(ax,colDesc,y,rows{q,3},'FontSize',10,'Color',[0 0 0], ...
                        'VerticalAlignment','top');
                case 'E'
                    et = eventTypes(rows{q,3});
                    if et.plain
                        patch(ax,colKey+[0 14 14 0],y+[3 3 15 15],[1 1 1],'EdgeColor',et.color);
                    else
                        patch(ax,colKey+[0 14 14 0],y+[3 3 15 15],et.color, ...
                            'EdgeColor',et.color,'FaceAlpha',0.6);
                    end
                    text(ax,colKey+22,y,rows{q,2},'FontWeight','bold','FontSize',10, ...
                        'Color',[0 0 0],'VerticalAlignment','top');
                    if et.plain
                        dsc = 'Transcript / notes only (no event, no colour)';
                    else
                        dsc = sprintf('Add / set "%s"  (%s outline)',et.name,styleName(et.style));
                    end
                    text(ax,colDesc,y,dsc,'FontSize',10,'Color',[0 0 0],'VerticalAlignment','top');
            end
            y = y + lineH;
        end
        uicontrol(sf,'Style','pushbutton','String','Close','Units','pixels', ...
            'Position',[(figW-110)/2 14 110 30],'Callback',@(~,~)delete(sf));

        function closeOnKey(src,e)
            if any(strcmp(e.Key,{'escape','return','space'})), delete(src); end
        end
    end

    % ---------------- helpers -----------------------------------------------
    function t = clampT(t), t = max(0,min(dur,t)); end

    function setStatus(msg,kind)
        if isempty(statusTxt) || ~ishghandle(statusTxt), return; end
        switch kind
            case 'warn', col = [0.75 0.35 0];
            case 'ok',   col = [0 0.45 0.10];
            otherwise,   col = [0.15 0.15 0.15];
        end
        set(statusTxt,'String',msg,'ForegroundColor',col);
    end

    function updateTitle()
        if isempty(fig) || ~ishghandle(fig), return; end
        if isempty(videoPath), set(fig,'Name',APPNAME); return; end
        [~,nm,ext] = fileparts(videoPath);
        star = ''; if dirty, star = '  * unsaved changes'; end
        set(fig,'Name',sprintf('%s - %s%s%s',APPNAME,nm,ext,star));
    end

    function s = fmtTime(t)
        if isempty(t) || isnan(t), s = '--:--.---'; return; end
        m = floor(t/60); s = sprintf('%02d:%06.3f',m,t-60*m);
    end

    function t = parseTime(s)
        % Accepts 14, 14.5, 14s, 1:05, 1:05.5, 0:01:05, 2m30, 2m30s, 2m.
        t = NaN; s = lower(strtrim(char(s))); if isempty(s), return; end
        s = strrep(s,' ','');
        s = regexprep(s,'(sec|secs|seconds|s)$','');
        if contains(s,'m')
            s = regexprep(s,'(min|mins|minutes|m)',':');
            if endsWith(s,':'), s = [s '0']; end
        end
        parts = strsplit(s,':'); v = str2double(parts);
        if any(isnan(v)) || numel(v) > 3 || any(v < 0), return; end
        t = 0; for q = 1:numel(v), t = t*60 + v(q); end
    end

    function s = styleName(ls)
        switch ls
            case '-',  s = 'solid';
            case '--', s = 'dashed';
            case ':',  s = 'dotted';
            otherwise, s = 'dash-dot';
        end
    end

    function deleteValid(h)
        if isempty(h), return; end
        h = h(ishghandle(h)); if ~isempty(h), delete(h); end
    end

    % uistack requires all objects to share a parent; restack each one.
    function stackEach(objs,where)
        objs = objs(ishghandle(objs));
        for o = reshape(objs,1,[])
            uistack(o,where);
        end
    end

    function x = toNum(v)
        if iscell(v), v = v{1}; end
        if isnumeric(v), x = double(v); else, x = str2double(string(v)); end
    end

    function s = toStr(v)
        if iscell(v), v = v{1}; end
        if isnumeric(v)
            if isnan(v), s = ''; else, s = num2str(v); end
        else
            s = char(string(v)); if strcmpi(s,'NaN')||strcmp(s,'<missing>'), s=''; end
        end
    end

    function onClose()
        try
            if ~confirmUnsaved('closing'), return; end
        catch
        end
        try, stop(playTimer); end %#ok<TRYNC>
        try, delete(playTimer); end %#ok<TRYNC>
        try, if ~isempty(mainPlayer), stop(mainPlayer); end; catch, end
        delete(fig);
    end
end
