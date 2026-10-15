function annotate_disfluencies_video(initialVideo)
% ANNOTATE_DISFLUENCIES_VIDEO  GUI for annotating stuttering-like disfluencies
% (SLD), typical disfluencies (TD) and speech errors from video.
%
%   annotate_disfluencies_video              opens the GUI (load a video with the button)
%   annotate_disfluencies_video(videoPath)   opens the GUI and loads the given video
%
% LAYOUT (top to bottom)
%   Annotation display + IPA buttons (left) | video
%   Spectrogram | Waveform | Transcript strip | Disfluency-Error strip | Notes | Overview
%
% WHAT YOU ANNOTATE
%   Highlight speech, double-click it, type what was said, e.g. "t-t-table".
%   It appears on the Transcript strip as pieces  t-  t-  table  (spread over
%   the highlight; drag their edges to fit the audio).
%   Double-click a highlight on the Notes strip to type notes for it.
%   Highlight a word + Enter = fluent event selected; highlight a part +
%   Enter = tick its disfluency / error type(s).
%   Right-click (two-finger click) a piece and tick what it is:
%     Fluent, Rep, Prol, Block, TD, Error > type, or Several...  (a piece can
%     have several). Disfluencies / errors connect to the fluent word they
%     overlap, or the next one after them. The same menu works on any
%     highlight on the spectrogram (e.g. a block over the whole episode).
%
% PLAYBACK
%   Drag to select: the selection stays put. Space plays it (pause keeps your
%   place), L loops it, 1-4 set the speed. A small bar on the selection offers
%   Play / Loop / speed / Create.
%
% Event times come from the AUDIO sample clock. Requires base MATLAB only.

% ----------------------------------------------------------------------------
% Configuration
% ----------------------------------------------------------------------------
APPNAME = 'Disfluency Annotator';
APPVER  = '4.1';

DISFL = struct( ...
    'abbr', {'rep',              'prol',               'bl',                 'td'}, ...
    'key',  {'r',                'p',                  'b',                  't'}, ...
    'name', {'repetition (SLD)', 'prolongation (SLD)', 'block (SLD)',        'typical disfluency'}, ...
    'color',{[0.85 0.15 0.15],   [0.10 0.60 0.20],     [0.15 0.35 0.90],     [0.95 0.78 0.05]});
ERRS = struct( ...
    'abbr', {'er-os',         'er-ow',        'er-is',          'er-iw',          'er-t',          'er-s'}, ...
    'name', {'omitted sound', 'omitted word', 'inserted sound', 'word insertion', 'transposition', 'substitution'});
ERR_KEY = 'e';   FLU_KEY = 'f';
FLUENT_COLOR = [0.78 0.78 0.78];
ERR_COLOR    = [0.85 0.00 0.85];
TRANS_COLOR  = [0.45 0.45 0.45];

% false (per spec): disfluency / error rows are saved WITHOUT a transcript, and
% a loaded file with transcripts on disfluency / error rows is flagged.
% (What was said still shows on screen from the transcript piece.)
% Set to true to save what was said ('t', 's' ...) on disfluency rows instead.
DISFL_HAS_TRANSCRIPT = false;

% IPA buttons, in this order (char codes keep the file encoding-safe):
%   row 1: θ ð ʃ ʒ x ŋ ɹ ʔ      row 2: ɛ æ ɪ ə ʊ ʌ ɑ ɔ ɚ ɝ
IPA_ROWS = {char([952 240 643 658 120 331 633 660]), ...
            char([603 230 618 601 650 652 593 596 602 605])};

ANN_COLOR  = [0 0 0];
SEL_COLOR  = [0.45 0.45 0.45];
OV_COLOR   = [0.20 0.45 0.95];
CUR_COLOR  = [0.90 0.50 0.00];
PLAY_COLOR = [0.00 0.55 0.55];
BG         = [0.94 0.94 0.94];
SPEEDS      = [1 0.75 0.5 0.25];
SPEED_NAMES = {'1x','0.75x','0.5x','0.25x'};

MINVIEW = 0.02;     % narrowest view (s)
MINSEL  = 0.005;    % shortest selection / event (s)
EDGE_PX = 8;        % grab tolerance for the resize handles (pixels)
MAXUNDO = 40;       % undo depth
TOL     = 1e-6;     % overlap tolerance (s) - touching edges are fine

% ----------------------------------------------------------------------------
% Shared state
% ----------------------------------------------------------------------------
initPath = '';
if nargin >= 1 && ~isempty(initialVideo), initPath = initialVideo; end

% kind: 'trans' (transcript text), 'fluent' / 'disfl' (parts), 'error'.
% link = uid of the fluent part a disfluency / error resolves to.
% wi   = (older files) word position; tok = uid of the transcript piece an event sits on.
EV0 = struct('start',{},'end',{},'kind',{},'cat',{},'transcript',{},'notes',{}, ...
    'uid',{},'link',{},'wi',{},'tok',{},'gStart',{},'gStop',{});

vid = []; videoPath = ''; audio = []; fs = 0; dur = 0; audioMax = 1; maxFreq = 5000;
viewStart = 0; viewEnd = 1; cursorTime = 0; selStart = NaN; selEnd = NaN;
events = EV0; nextUid = 1;
currentEvent = []; dirty = false; undoStack = {}; listMap = [];

pendingNew = [];             % struct(kind,a,b) while typing a new transcript
batchUids = [];              % events added together with M share the typed text
laneOf = []; nLaneP = 1; nLaneE = 1;   % row (lane) of each item on the Disfl/Error strip
lastKeyCommit = 0;           % time (s) of last Enter/Esc handled by a text box
lastBoxKey = '';             % last key typed in a text box (Enter commits)

isEditing = false; typingBox = []; transBox = []; notesBox = []; deBox = [];
isPlaying = false; playStartTime = 0; playEndTime = 0; playSpeed = 1; mainPlayer = [];
loopOn = true; playLoopA = NaN; loopPass = 0;
dragAxes = []; dragMode = ''; dragStartT = 0; dragY = NaN; dragPrevCur = [];
downPix = [0 0]; didDrag = false; dragPanel = ''; moveIdx = []; moveOrig = [0 0];
moveKids = []; moveKidsOrig = zeros(2,0); movePushed = false;
resizeUndoPushed = false; frameStride = 1;

imgVideo = []; placeholderTxt = []; hSpecImg = []; hWave = [];
evGfx = gobjects(0); selGfx = gobjects(0); curGfx = gobjects(0); playGfx = gobjects(0);
ovEnvGfx = gobjects(0); ovDyn = gobjects(0); ovCur = gobjects(0);
ovMode = ''; ovGrab = 0; edtFrom = []; edtTo = []; btnGoTF = [];

fig = []; statusTxt = []; infoTxt = []; viewTxt = []; lstEvents = []; annotDisp = [];
popSpeed = []; chkLoop = []; videoDependent = [];
selBar = []; btnBarPlay = []; chkBarLoop = []; popBarSpeed = []; SELBAR_W = 188;
vzoom = []; lastFull = [];               % video zoom region [x0 x1 y0 y1]; last full frame
fcFrames = []; fcTimes = []; fcRange = []; fcZoom = []; fcStride = 1; fcFps = 30; fcLast = 0;
vpanStart = [0 0]; vpanZoom = [];
groupU = 0;                          % uid of the selected fluent word (Enter workflow)
expandedU = []; listKind = {};       % events table: open fluent rows, row types
unresolvedOpen = false;              % events table: "(unresolved)" section open?
popFig = []; vidPosNorm = []; dockMsg = []; btnPop = []; vidBtns = gobjects(1,3);   % pop-out video window
CACHE_MAX_S  = 60;     % longest highlight decoded into memory for smooth playback (s)
CACHE_MAX_MB = 1200;   % memory cap for that (MB)
% Smooth playback copy: on load, a video that is taller than SMOOTH_HEIGHT or
% not .mp4 (e.g. iPhone .MOV / HEVC) is converted once with ffmpeg to H.264,
% constant frame rate, no audio, saved next to it as <name>_smooth.mp4 and
% reused after that. Audio (and so every annotation time) still comes from the
% ORIGINAL file. Needs ffmpeg (Mac: brew install ffmpeg). false = never.
SMOOTH_ENABLE = true;
SMOOTH_HEIGHT = 720;
SMOOTH_SUFFIX = '_smooth';
cmenu = []; miFlu = []; miDis = []; miErrRoot = []; miErr = []; miSev = []; miSplit = [];
miEdit = []; miPlay = []; miDel = []; resizeNbr = [];
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

% annotation display (left of video): full text of the selected item
pnlAnn = uipanel(fig,'Title','Selected annotation','Units','normalized', ...
    'Position',[0.04 0.668 0.225 0.307],'BackgroundColor',BG,'FontWeight','bold');
annotDisp = uicontrol(pnlAnn,'Style','edit','Max',2,'Min',0,'Enable','inactive', ...
    'Units','normalized','Position',[0.03 0.03 0.94 0.94],'HorizontalAlignment','left', ...
    'FontSize',11,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
    'String',{'Nothing selected.'});

% IPA buttons: click to insert into the text box you are typing in
pnlIPA = uipanel(fig,'Units','normalized','Position',[0.04 0.600 0.225 0.064], ...
    'BackgroundColor',BG,'BorderType','none');
for ipaRow = 1:numel(IPA_ROWS)
    ipaSyms = IPA_ROWS{ipaRow};
    for ipaI = 1:numel(ipaSyms)
        uicontrol(pnlIPA,'Style','pushbutton','String',ipaSyms(ipaI),'Units','normalized', ...
            'Position',[(ipaI-1)*0.1 (numel(IPA_ROWS)-ipaRow)*0.5+0.03 0.096 0.44], ...
            'FontSize',12,'TooltipString','Insert into the text you are typing', ...
            'Callback',@(src,~)cb(@()insertSymbol(src)));
    end
end

axVideo = axes('Parent',fig,'Units','normalized','Position',[0.275 0.600 0.465 0.375], ...
    'Color',[0 0 0],'SortMethod','childorder');
axis(axVideo,[0 1 0 1]); axis(axVideo,'off');
placeholderTxt = text(axVideo,0.5,0.5,'Click "Load video" to begin', ...
    'HorizontalAlignment','center','FontSize',14,'Color',[0.6 0.6 0.6]);
% video zoom: scroll on the picture, drag to move, double-click / Fit to reset
set(axVideo,'PickableParts','all','ButtonDownFcn',@(~,~)cb(@onVideoDown));
btnPop = uicontrol(fig,'Style','pushbutton','String','Pop out video','Units','normalized', ...
    'FontSize',9,'Position',[0.588 0.948 0.070 0.025], ...
    'TooltipString','Open the video in its own window (resize it for a bigger, sharper picture)', ...
    'Callback',@(src,~)cb(@()barDo(src,@togglePopVideo)));
vidBtns(1) = uicontrol(fig,'Style','pushbutton','String','+','Units','normalized','FontSize',12, ...
    'Position',[0.660 0.948 0.022 0.025],'TooltipString','Zoom the video in (or scroll on it)', ...
    'Callback',@(src,~)cb(@()barDo(src,@()videoZoomBy(1.5))));
vidBtns(2) = uicontrol(fig,'Style','pushbutton','String','-','Units','normalized','FontSize',12, ...
    'Position',[0.684 0.948 0.022 0.025],'TooltipString','Zoom the video out', ...
    'Callback',@(src,~)cb(@()barDo(src,@()videoZoomBy(1/1.5))));
vidBtns(3) = uicontrol(fig,'Style','pushbutton','String','Fit','Units','normalized','FontSize',9, ...
    'Position',[0.708 0.948 0.030 0.025],'TooltipString','Show the whole picture (or double-click it)', ...
    'Callback',@(src,~)cb(@()barDo(src,@()resetVideoZoom())));

axSpec     = axes('Parent',fig,'Units','normalized','Position',[0.07 0.465 0.67 0.125]);
axWave     = axes('Parent',fig,'Units','normalized','Position',[0.07 0.397 0.67 0.062]);
axTrans    = axes('Parent',fig,'Units','normalized','Position',[0.07 0.345 0.67 0.046]);
axDE       = axes('Parent',fig,'Units','normalized','Position',[0.07 0.267 0.67 0.072]);
axNotes    = axes('Parent',fig,'Units','normalized','Position',[0.07 0.222 0.67 0.039]);
axOverview = axes('Parent',fig,'Units','normalized','Position',[0.07 0.138 0.67 0.036]);

for axInit = [axVideo axSpec axWave axTrans axDE axNotes axOverview]
    try, disableDefaultInteractivity(axInit); catch, end
    try, axInit.Toolbar.Visible = 'off'; catch, end
end
for axInit = [axSpec axWave axTrans axDE axNotes axOverview]
    hold(axInit,'on'); box(axInit,'on');
    set(axInit,'XTick',[],'YTick',[],'Color',[1 1 1], ...
        'XColor',[0.15 0.15 0.15],'YColor',[0.15 0.15 0.15]);
end
set(axSpec,'YDir','normal','YTickMode','auto'); ylabel(axSpec,'Freq (Hz)');
ylabel(axWave,'Amp'); ylim(axWave,[-1 1]);
ylabel(axTrans,'Transcript','Rotation',0,'HorizontalAlignment','right'); ylim(axTrans,[0 1]);
ylabel(axDE,{'Disfl /','Error'},'Rotation',0,'HorizontalAlignment','right'); ylim(axDE,[0 1]);
ylabel(axNotes,'Notes','Rotation',0,'HorizontalAlignment','right');      ylim(axNotes,[0 1]);
ylabel(axOverview,'Overview','Rotation',0,'HorizontalAlignment','right'); ylim(axOverview,[0 1]);

set(axSpec,    'ButtonDownFcn',@(~,~)cb(@()onAxDown(axSpec)));
set(axWave,    'ButtonDownFcn',@(~,~)cb(@()onAxDown(axWave)));
set(axTrans,   'ButtonDownFcn',@(~,~)cb(@()onAxDown(axTrans,'transcript')));
set(axDE,      'ButtonDownFcn',@(~,~)cb(@()onAxDown(axDE,'disferr')));
set(axNotes,   'ButtonDownFcn',@(~,~)cb(@()onAxDown(axNotes,'notes')));
set(axOverview,'ButtonDownFcn',@(~,~)cb(@onOverviewDown));

% Inline text boxes shown over an item while you type. Enter saves, Esc cancels.
transBox = uicontrol(fig,'Style','edit','Max',1,'Min',0,'Units','pixels', ...
    'HorizontalAlignment','left','FontSize',10,'Visible','off','Enable','inactive', ...
    'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Tag','transcript', ...
    'TooltipString','Type the transcript - Enter saves, Esc cancels', ...
    'ButtonDownFcn',@(~,~)cb(@()beginTyping('text')), ...
    'KeyPressFcn',@(~,e)cb(@()boxKey(e)), ...
    'Callback',@(src,~)cb(@()commitBox(src,false)));
deBox = uicontrol(fig,'Style','edit','Max',1,'Min',0,'Units','pixels', ...
    'HorizontalAlignment','left','FontSize',10,'Visible','off','Enable','inactive', ...
    'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Tag','detext', ...
    'TooltipString','Type what was said - Enter saves, Esc cancels', ...
    'ButtonDownFcn',@(~,~)cb(@()beginTyping('text')), ...
    'KeyPressFcn',@(~,e)cb(@()boxKey(e)), ...
    'Callback',@(src,~)cb(@()commitBox(src,false)));
notesBox = uicontrol(fig,'Style','edit','Max',1,'Min',0,'Units','pixels', ...
    'HorizontalAlignment','left','FontSize',10,'Visible','off','Enable','inactive', ...
    'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Tag','notes', ...
    'TooltipString','Type notes - Enter saves, Esc cancels', ...
    'ButtonDownFcn',@(~,~)cb(@()beginTyping('notes')), ...
    'KeyPressFcn',@(~,e)cb(@()boxKey(e)), ...
    'Callback',@(src,~)cb(@()commitBox(src,false)));

% Selection bar: sits on the selection - replay it, loop it, change speed, or create
selBar = uipanel(fig,'Units','pixels','Position',[0 0 344 30],'BackgroundColor',[1 1 1], ...
    'BorderType','line','HighlightColor',[0.3 0.3 0.3],'Visible','off');
btnBarPlay = uicontrol(selBar,'Style','pushbutton','String','Play','Units','pixels', ...
    'Position',[3 3 52 22],'TooltipString','Play / pause the selection (Space)', ...
    'Callback',@(src,~)cb(@()barDo(src,@togglePlay)));
chkBarLoop = uicontrol(selBar,'Style','checkbox','String','Loop','Value',1,'Units','pixels', ...
    'Position',[59 3 52 22],'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
    'TooltipString','Replay the selection until you pause (L)', ...
    'Callback',@(src,~)cb(@()barDo(src,@()setLoop(get(src,'Value')==1))));
popBarSpeed = uicontrol(selBar,'Style','popupmenu','String',SPEED_NAMES,'Value',1, ...
    'Units','pixels','Position',[113 3 70 22],'TooltipString','Playback speed (keys 1-4)', ...
    'Callback',@(src,~)cb(@()barDo(src,@()setSpeedIdx(get(src,'Value')))));
SELBAR_W = 188;

% right-click (two-finger click) menu for a transcript piece or a highlight
cmenu = uicontextmenu(fig);
try, set(cmenu,'ContextMenuOpeningFcn',@(src,~)cb(@()onMenuOpen(src)));
catch, set(cmenu,'Callback',@(src,~)cb(@()onMenuOpen(src))); end
miFlu = uimenu(cmenu,'Label','Fluent','Callback',@(~,~)cb(@()menuToggle('fluent','')));
miDis = gobjects(1,numel(DISFL));
for miInit = 1:numel(DISFL)
    miDis(miInit) = uimenu(cmenu,'Label',sprintf('%s  (%s)',DISFL(miInit).name,DISFL(miInit).abbr), ...
        'Callback',@(~,~)cb(@()menuToggle('disfl',DISFL(miInit).abbr)));
end
miErrRoot = uimenu(cmenu,'Label','Error');
miErr = gobjects(1,numel(ERRS));
for miInit = 1:numel(ERRS)
    miErr(miInit) = uimenu(miErrRoot,'Label',sprintf('%s  -  %s',ERRS(miInit).abbr,ERRS(miInit).name), ...
        'Callback',@(~,~)cb(@()menuToggle('error',ERRS(miInit).abbr)));
end
miSev   = uimenu(cmenu,'Label','Several...','Callback',@(~,~)cb(@onMultiKey));
miSplit = uimenu(cmenu,'Label','Split into words','Separator','on','Callback',@(~,~)cb(@menuSplit));
miEdit  = uimenu(cmenu,'Label','Edit text','Callback',@(~,~)cb(@editCurrentEvent));
miPlay  = uimenu(cmenu,'Label','Play','Callback',@(~,~)cb(@()startPlayFresh()));
miDel   = uimenu(cmenu,'Label','Delete','Callback',@(~,~)cb(@deleteCurrentEvent));
for axInit = [axSpec axWave axTrans axDE axNotes]
    set(axInit,'UIContextMenu',cmenu);
end

splashStep(0.55,'Creating controls...');

% ---- right column: events panel ----------------------------------------------
pnlEv = uipanel(fig,'Title','Fluent events (click to open)','Units','normalized', ...
    'Position',[0.765 0.30 0.22 0.675],'BackgroundColor',BG,'FontWeight','bold');
lstEvents = uitable(pnlEv,'Units','normalized','Position',[0.04 0.40 0.92 0.58], ...
    'ColumnName',{'','Start','End','Fluent event','Notes'}, ...
    'ColumnWidth',{16 66 66 74 70},'ColumnEditable',false(1,5),'RowName',[], ...
    'Data',cell(0,5),'FontSize',9,'ForegroundColor',[0 0 0], ...
    'BackgroundColor',[1 1 1],'RowStriping','on', ...
    'CellSelectionCallback',@(src,e)cb(@()onListSelect(src,e)));
btnPlayEv = mkBtn(pnlEv,'Play event',  [0.04 0.320 0.45 0.065],@playCurrentEvent, ...
    'Play the selected item');
btnEditEv = mkBtn(pnlEv,'Edit text',   [0.51 0.320 0.45 0.065],@editCurrentEvent, ...
    'Type the text of the selected item (Enter)');
btnTypeEv = mkBtn(pnlEv,'Change type', [0.04 0.245 0.45 0.065],@changeTypeDialog, ...
    'Change the type of the selected part / error');
btnDelEv  = mkBtn(pnlEv,'Delete event',[0.51 0.245 0.45 0.065],@deleteCurrentEvent, ...
    'Delete the selected item (Ctrl+Z to undo)');
btnAddTxt = mkBtn(pnlEv,'Transcript',[0.04 0.170 0.45 0.065],@()beginNewText('trans'), ...
    'Type the transcript for the selection (or double-click the highlight)');
btnAddFlu = mkBtn(pnlEv,'Fluent part (F)',[0.51 0.170 0.45 0.065],@()onPartKey('flu'), ...
    'Mark the selection as the fluent word (inside a transcript)');
uicontrol(pnlEv,'Style','text','String','Playback speed:','Units','normalized', ...
    'Position',[0.04 0.100 0.45 0.045],'HorizontalAlignment','left','BackgroundColor',BG);
popSpeed = uicontrol(pnlEv,'Style','popupmenu','Units','normalized', ...
    'String',SPEED_NAMES,'Value',1, ...
    'Position',[0.51 0.105 0.45 0.050],'Callback',@(src,~)cb(@()onSpeed(src)));
infoTxt = uicontrol(pnlEv,'Style','text','String','No selection','Units','normalized', ...
    'Position',[0.04 0.010 0.92 0.085],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontSize',8.5);

pnlTips = uipanel(fig,'Title','Quick tips','Units','normalized', ...
    'Position',[0.765 0.135 0.22 0.155],'BackgroundColor',BG,'FontWeight','bold');
uicontrol(pnlTips,'Style','text','Units','normalized','Position',[0.03 0.02 0.94 0.96], ...
    'HorizontalAlignment','left','BackgroundColor',BG,'FontSize',8.5,'String', { ...
    'New speech: highlight + double-click = transcript.', ...
    'Notes: double-click highlight on Notes strip.', ...
    'Word: highlight + Enter = fluent event selected.', ...
    'Part: highlight + Enter = tick type(s), Enter.', ...
    'W = whole episode. Esc = done. N = notes.', ...
    'Space plays (L loop, 1-4 speed). F1 = help'});

% ---- navigation row ------------------------------------------------------------
uicontrol(fig,'Style','text','String','Navigate:','Units','normalized', ...
    'Position',[0.04 0.083 0.055 0.035],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontWeight','bold');
btnZin   = mkBtn(fig,'Zoom in',        [0.095 0.087 0.065 0.040],@()zoomAbout(cursorCenter(),0.5), ...
    'Zoom in around the cursor (Ctrl+I)');
btnZout  = mkBtn(fig,'Zoom out',       [0.165 0.087 0.065 0.040],@()zoomAbout(cursorCenter(),2), ...
    'Zoom out around the cursor (Ctrl+O)');
btnZsel  = mkBtn(fig,'Fit selection',  [0.235 0.087 0.080 0.040],@zoomToSelection, ...
    'Zoom to the current selection (Ctrl+N)');
btnZall  = mkBtn(fig,'Full recording', [0.320 0.087 0.085 0.040],@fullView, ...
    'Show the whole recording (Ctrl+A)');
chkLoop  = uicontrol(fig,'Style','checkbox','String','Loop (L)','Value',1,'Units','normalized', ...
    'Position',[0.412 0.087 0.065 0.040],'BackgroundColor',BG,'ForegroundColor',[0 0 0], ...
    'TooltipString','Replay the selection over and over until you pause', ...
    'Callback',@(src,~)cb(@()setLoop(get(src,'Value')==1)));
viewTxt = uicontrol(fig,'Style','text','String','','Units','normalized', ...
    'Position',[0.07 0.176 0.36 0.024],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontSize',9);
uicontrol(fig,'Style','text','String','Show from','Units','normalized', ...
    'Position',[0.432 0.176 0.055 0.024],'HorizontalAlignment','right', ...
    'BackgroundColor',BG,'FontSize',9);
edtFrom = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.490 0.174 0.085 0.028], ...
    'FontSize',9,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Enable','inactive', ...
    'Tag','tfrom','HorizontalAlignment','center', ...
    'TooltipString','Start time (seconds or m:ss) - click, type, Enter', ...
    'ButtonDownFcn',@(src,~)cb(@()beginFieldEdit(src)), ...
    'KeyPressFcn',@(~,e)cb(@()boxKey(e)),'Callback',@(src,~)cb(@()commitBox(src)));
uicontrol(fig,'Style','text','String','to','Units','normalized', ...
    'Position',[0.577 0.176 0.018 0.024],'HorizontalAlignment','center', ...
    'BackgroundColor',BG,'FontSize',9);
edtTo = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.597 0.174 0.085 0.028], ...
    'FontSize',9,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Enable','inactive', ...
    'Tag','tto','HorizontalAlignment','center', ...
    'TooltipString','End time (seconds or m:ss) - click, type, Enter', ...
    'ButtonDownFcn',@(src,~)cb(@()beginFieldEdit(src)), ...
    'KeyPressFcn',@(~,e)cb(@()boxKey(e)),'Callback',@(src,~)cb(@()commitBox(src)));
btnGoTF = mkBtn(fig,'Go',[0.686 0.173 0.054 0.030],@applyTimeFields, ...
    'Show the time frame typed in the From / To boxes');

% ---- file row --------------------------------------------------------------------
mkBtn(fig,'Load video',            [0.040 0.035 0.100 0.042],@uiLoadVideo, ...
    'Open a video file');
btnLoadAnnot = mkBtn(fig,'Load annotations...',[0.145 0.035 0.120 0.042],@loadAnnotations, ...
    'Load an annotation table (checks it first); merge or replace');
btnSave  = mkBtn(fig,'Save annotations',[0.270 0.035 0.120 0.042],@saveAnnotations, ...
    'Save all annotations to an Excel file (Ctrl+S)');
btnUndo  = mkBtn(fig,'Undo',            [0.395 0.035 0.070 0.042],@undo, ...
    'Undo the last change (Ctrl+Z)');
mkBtn(fig,'Help',                  [0.470 0.035 0.070 0.042],@showHelp, ...
    'Full instructions (F1)');
mkBtn(fig,'Shortcuts',             [0.545 0.035 0.080 0.042],@showShortcuts, ...
    'Keyboard and mouse shortcuts');

% colour legend
legNames = [{sprintf('%s fluent',upper(FLU_KEY))}, ...
    arrayfun(@(d)sprintf('%s %s',upper(d.key),d.abbr),DISFL,'UniformOutput',false), ...
    {sprintf('%s error',upper(ERR_KEY))}];
legCols  = [{FLUENT_COLOR}, {DISFL.color}, {ERR_COLOR}];
legX = 0.630;
for liInit = 1:numel(legNames)
    uicontrol(fig,'Style','text','String','','Units','normalized', ...
        'Position',[legX 0.047 0.009 0.018],'BackgroundColor',legCols{liInit});
    uicontrol(fig,'Style','text','Units','normalized','Position',[legX+0.011 0.038 0.041 0.032], ...
        'String',legNames{liInit},'HorizontalAlignment','left','BackgroundColor',BG,'FontSize',8);
    legX = legX + 0.058;
end

statusTxt = uicontrol(fig,'Style','text','String','No video loaded - click "Load video" to begin.', ...
    'Units','normalized','Position',[0.04 0.004 0.94 0.026],'HorizontalAlignment','left', ...
    'BackgroundColor',BG,'FontSize',9,'ForegroundColor',[0.2 0.2 0.2]);

set(findall(fig,'Style','text'),'ForegroundColor',[0 0 0]);
set(findall(fig,'Type','uipanel'),'ForegroundColor',[0 0 0]);

videoDependent = [edtFrom edtTo btnGoTF lstEvents btnPlayEv btnEditEv btnTypeEv btnDelEv ...
    btnAddTxt btnAddFlu popSpeed btnZin btnZout btnZsel btnZall chkLoop ...
    btnLoadAnnot btnSave btnUndo btnPop];
set(videoDependent,'Enable','off');

splashStep(0.80,'Setting up playback engine...');
playTimer = timer('ExecutionMode','fixedRate','Period',0.02,'BusyMode','drop', ...
    'TimerFcn',@(~,~)cb(@playTick));

splashStep(1.00,'Ready');
closeSplash();
set(fig,'Visible','on'); drawnow;

showWelcome();
if ~isempty(initPath), cb(@()loadVideo(initPath)); end

% ============================================================================
%                               NESTED FUNCTIONS
% ============================================================================

    % ======================== annotation model (v5) ===========================

    function showWelcome()
        try
            if ispref('DisfluencyAnnotator','hideWelcome') && ...
                    getpref('DisfluencyAnnotator','hideWelcome')
                return;
            end
        catch
        end
        msg = sprintf([ ...
            'How it works (example "the t-t-table is on the ground"):\n\n' ...
            '1. Highlight the speech, DOUBLE-CLICK it, type the sentence, Enter.\n' ...
            '   It appears on the Transcript strip as pieces: the  t-  t-  table ...\n' ...
            '   (Double-click a highlight on the Notes strip to type notes.)\n' ...
            '2. Highlight "table" (or click its piece) and press Enter:\n' ...
            '   FLUENT EVENT SELECTED.\n' ...
            '3. Highlight "t-t" (or one "t-", or part of a word) and press Enter.\n' ...
            '   Tick the type(s): rep / prol / bl / td and / or errors (keys\n' ...
            '   R P B T, 1-6). Enter to finish that part. Repeat for other parts.\n' ...
            '   W = the whole episode (t-t-table) at once.\n' ...
            '4. Esc when you are done with that word (or just double-click the\n' ...
            '   next highlight to type its transcript).\n\n' ...
            'Space plays the highlight (L loop, 1-4 speed). N = notes. Ctrl+S saves.\n' ...
            'F1 or Help = full instructions.']);
        c = questdlg(msg,['Welcome to ' APPNAME],'Get started','Don''t show again','Get started');
        if strcmp(c,'Don''t show again')
            try, setpref('DisfluencyAnnotator','hideWelcome',true); catch, end
        end
    end

    function redrawEvents()
        deleteValid(evGfx); evGfx = gobjects(0);
        if isempty(vid), return; end
        computeLanes();
        vw = viewEnd - viewStart;
        yl = get(axSpec,'YLim'); ys = diff(yl);
        selU = 0;
        if ~isempty(currentEvent) && currentEvent <= numel(events), selU = groupUid(currentEvent); end
        for k = 1:numel(events)
            ev = events(k);
            if ev.end < viewStart || ev.start > viewEnd, continue; end
            col = eventColor(k);
            isSel = isequal(currentEvent,k);
            if isSel, fa = 0.40; lw = 2; else, fa = 0.22; lw = 1; end
            xs = [ev.start ev.end ev.end ev.start];
            tx0 = max(ev.start,viewStart) + 0.003*vw;
            g = gobjects(0);
            switch ev.kind
                case 'trans'
                    % a transcript piece, striped with EVERY annotation tied to it:
                    % its own marks, and for a fluent word also the disfluencies /
                    % errors that resolve to it
                    pc = pieceColors(k);
                    if isempty(pc)
                        g(end+1) = patch(axTrans,'XData',xs,'YData',[0.04 0.04 0.96 0.96], ...
                            'FaceColor',[0.96 0.96 0.96],'EdgeColor',[0.35 0.35 0.35],'LineWidth',lw); %#ok<AGROW>
                    else
                        nb = size(pc,1); hb = 0.92/nb;
                        for q = 1:nb
                            y0 = 0.96 - q*hb;
                            g(end+1) = patch(axTrans,'XData',xs,'YData',[y0 y0 y0+hb y0+hb], ...
                                'FaceColor',0.45 + 0.55*pc(q,:),'EdgeColor','none'); %#ok<AGROW>
                        end
                        g(end+1) = patch(axTrans,'XData',xs,'YData',[0.04 0.04 0.96 0.96], ...
                            'FaceColor','none','EdgeColor',[0.35 0.35 0.35],'LineWidth',lw); %#ok<AGROW>
                    end
                    isFlu = any(strcmp({events(attachedTo(k)).kind},'fluent'));
                    inSel = selU > 0 && any([events(attachedTo(k)).uid] == selU | ...
                        [events(attachedTo(k)).link] == selU);
                    if inSel, tc = [0.85 0 0]; elseif isFlu, tc = [0.10 0.30 0.85]; else, tc = ANN_COLOR; end
                    if isFlu || inSel, fw = 'bold'; else, fw = 'normal'; end
                    g(end+1) = text(axTrans,tx0,0.5,fitText(ev.transcript,axTrans,ev,10),'Clipping','on', ...
                        'Interpreter','none','FontSize',10,'FontWeight',fw,'Color',tc, ...
                        'VerticalAlignment','middle'); %#ok<AGROW>
                case 'fluent'
                    % shown by its transcript piece (blue word) - no box of its own
                otherwise
                    if strcmp(ev.kind,'error')
                        g(end+1) = patch(axSpec,'XData',xs,'YData',[yl(2)-0.12*ys yl(2)-0.12*ys yl(2) yl(2)], ...
                            'FaceColor',col,'FaceAlpha',0.55,'EdgeColor',col,'LineWidth',lw); %#ok<AGROW>
                        g(end+1) = patch(axWave,'XData',xs,'YData',[-1 -1 -0.76 -0.76], ...
                            'FaceColor',col,'FaceAlpha',0.55,'EdgeColor',col,'LineWidth',lw); %#ok<AGROW>
                        faD = 0.35;
                    else
                        if strcmp(ev.kind,'fluent'), faS = 0.30; else, faS = fa; end
                        g(end+1) = patch(axSpec,'XData',xs,'YData',[yl(1) yl(1) yl(2) yl(2)], ...
                            'FaceColor',col,'FaceAlpha',faS,'EdgeColor',col*0.8,'LineWidth',lw); %#ok<AGROW>
                        g(end+1) = patch(axWave,'XData',xs,'YData',[-1 -1 1 1], ...
                            'FaceColor',col,'FaceAlpha',faS,'EdgeColor',col*0.8,'LineWidth',lw); %#ok<AGROW>
                        faD = max(faS,0.35);
                    end
                    bnd = laneBand(k);
                    g(end+1) = patch(axDE,'XData',xs,'YData',[bnd(1) bnd(1) bnd(2) bnd(2)], ...
                        'FaceColor',col,'FaceAlpha',faD,'EdgeColor',col*0.8,'LineWidth',lw); %#ok<AGROW>
                    if diff(bnd) < 0.11, fsz = 7; else, fsz = 8; end
                    g(end+1) = text(axDE,tx0,mean(bnd),fitText(deLabel(k),axDE,ev,fsz),'Clipping','on','Interpreter','none', ...
                        'FontSize',fsz,'FontWeight','bold','Color',ANN_COLOR, ...
                        'VerticalAlignment','middle'); %#ok<AGROW>
            end
            if ~isempty(ev.notes)
                switch ev.kind
                    case 'trans', ty = 0.78;
                    case 'error', ty = 0.22;
                    otherwise,    ty = 0.5;
                end
                g(end+1) = patch(axNotes,'XData',xs,'YData',[ty-0.2 ty-0.2 ty+0.2 ty+0.2], ...
                    'FaceColor','none','EdgeColor',col,'LineWidth',lw); %#ok<AGROW>
                g(end+1) = text(axNotes,tx0,ty,fitText(ev.notes,axNotes,ev,8),'Clipping','on','Interpreter','none', ...
                    'FontSize',8,'Color',ANN_COLOR,'VerticalAlignment','middle'); %#ok<AGROW>
            end
            evGfx = [evGfx g]; %#ok<AGROW>
        end

        % bold black boxes around the rest of the selected word's group
        for k = focusSet()
            ev = events(k);
            if ev.end < viewStart || ev.start > viewEnd, continue; end
            xs = [ev.start ev.end ev.end ev.start];
            if strcmp(ev.kind,'trans')
                evGfx(end+1) = patch(axTrans,'XData',xs,'YData',[0 0 1 1], ...
                    'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',3); %#ok<AGROW>
            elseif ~strcmp(ev.kind,'fluent')
                bnd = laneBand(k);
                evGfx(end+1) = patch(axDE,'XData',xs,'YData',[bnd(1) bnd(1) bnd(2) bnd(2)], ...
                    'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',2.5); %#ok<AGROW>
            end
        end

        % group mode: label the selected fluent word
        if groupU > 0
            f = findUid(groupU);
            if ~isempty(f) && events(f).end > viewStart && events(f).start < viewEnd
                xs = [events(f).start events(f).end events(f).end events(f).start];
                for hx = [axSpec axWave]
                    yl2 = get(hx,'YLim');
                    evGfx(end+1) = patch(hx,'XData',xs,'YData',[yl2(1) yl2(1) yl2(2) yl2(2)], ...
                        'FaceColor',FLUENT_COLOR,'FaceAlpha',0.25,'EdgeColor',[0 0 0],'LineWidth',2.5); %#ok<AGROW>
                end
                xm = (max(events(f).start,viewStart) + min(events(f).end,viewEnd))/2;
                evGfx(end+1) = text(axSpec,xm,yl(1)+0.5*ys,'fluent event selected', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10, ...
                    'FontWeight','bold','Color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',2, ...
                    'Clipping','on','Interpreter','none'); %#ok<AGROW>
                evGfx(end+1) = text(axWave,xm,0,'highlight a part + Enter  |  Esc = done', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8, ...
                    'Color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',1, ...
                    'Clipping','on','Interpreter','none'); %#ok<AGROW>
            end
        end

        if ~isempty(evGfx), set(evGfx,'PickableParts','none','HitTest','off'); end
        drawOverview();
        positionBoxes();
        updateAnnotDisplay();
    end

    function s = transTex(sk)
        % Transcript text with annotated words highlighted: blue = annotated
        % word, red = the word whose group is selected.
        toks = regexp(events(sk).transcript,'\S+','match');
        if isempty(toks), s = ''; return; end
        mark = zeros(1,numel(toks));
        selU = 0;
        if ~isempty(currentEvent) && currentEvent <= numel(events), selU = groupUid(currentEvent); end
        for f = fluentsInSeg(sk)
            wi = events(f).wi;
            if wi >= 1 && wi <= numel(toks)
                mark(wi) = max(mark(wi), 1 + (events(f).uid == selU));
            end
        end
        parts = cell(1,numel(toks));
        for i = 1:numel(toks)
            t = regexprep(toks{i},'([\\{}_^])','\\$1');
            switch mark(i)
                case 1, parts{i} = ['{\bf\color[rgb]{0.10 0.30 0.85}' t '}'];
                case 2, parts{i} = ['{\bf\color[rgb]{0.85 0.00 0.00}' t '}'];
                otherwise, parts{i} = t;
            end
        end
        s = strjoin(parts,' ');
    end

    function idx = focusSet()
        % Outlined in bold black: the rest of the selected item's word group
        % (fluent word + its disfluencies / errors) and their transcript pieces.
        idx = [];
        if isempty(events), return; end
        if groupU > 0 && ~isempty(findUid(groupU))
            u = groupU;
        else
            if isempty(currentEvent) || currentEvent > numel(events), return; end
            u = groupUid(currentEvent);
        end
        if u <= 0, return; end
        m = [findUid(u) linkedAll(u)];
        tk = [events(m).tok]; tk = tk(tk > 0);
        for t = unique(tk)
            j = findUid(t); if ~isempty(j), m(end+1) = j; end %#ok<AGROW>
        end
        idx = reshape(setdiff(m,currentEvent),1,[]);
    end

    function w = normWord(s)
        w = lower(strtrim(char(s)));
        w = regexprep(w,'^[\.,;:!?"''()\[\]-]+|[\.,;:!?"''()\[\]-]+$','');
    end

    function col = eventColor(k)
        switch events(k).kind
            case 'fluent', col = FLUENT_COLOR;
            case 'error',  col = ERR_COLOR;
            case 'trans',  col = TRANS_COLOR;
            otherwise
                q = find(strcmp({DISFL.abbr},events(k).cat),1);
                if isempty(q), col = [0.5 0.5 0.5]; else, col = DISFL(q).color; end
        end
    end

    function updateSelectionGraphics()
        deleteValid(selGfx); selGfx = gobjects(0);
        updateInfo();
        if ~hasSel() || isempty(vid), updateSelBar(); return; end
        for hx = [axSpec axWave axTrans axDE axNotes]
            yl = get(hx,'YLim');
            p = patch(hx,'XData',[selStart selEnd selEnd selStart], ...
                'YData',[yl(1) yl(1) yl(2) yl(2)],'FaceColor',SEL_COLOR, ...
                'FaceAlpha',0.12,'EdgeColor','none');
            selGfx = [selGfx p drawHandles(hx,selStart,selEnd)]; %#ok<AGROW>
        end
        set(selGfx,'PickableParts','none','HitTest','off'); stackEach(selGfx,'top');
        updateSelBar();
    end

    function updateSelBar()
        % Floating Play / Loop / speed / Create bar at the top of the selection.
        if isempty(selBar) || ~ishghandle(selBar), return; end
        % stays visible while you type a syllable's text: clicking another
        % button saves the text and adds the next event
        if isempty(vid) || ~hasSel() || ~isempty(pendingNew)
            set(selBar,'Visible','off'); return;
        end
        axpix = getpixelposition(axSpec,true); xl = get(axSpec,'XLim');
        w = SELBAR_W; h = 30;
        mid = (max(selStart,xl(1)) + min(selEnd,xl(2)))/2;
        f = (mid - xl(1))/max(diff(xl),eps);
        px = axpix(1) + f*axpix(3) - w/2;
        px = max(min(px, axpix(1)+axpix(3)-w), axpix(1));
        py = axpix(2) + axpix(4) - h - 2;
        set(selBar,'Position',[px py w h],'Visible','on');
        syncPlayBtn();
    end

    function barDo(src,fn)
        dropFocus(src);
        if isEditing, commitEdit(); end
        fn();
    end

    function dropFocus(src)
        % Toggling Enable hands keyboard focus back to the figure (Space, keys).
        try
            set(src,'Enable','off'); drawnow; set(src,'Enable','on');
        catch
        end
    end

    function barAdd(src)
        v = get(src,'Value'); set(src,'Value',1);
        switch v
            case 2, beginNewText('trans');
            case 3, onPartKey('flu');
            case {4,5,6,7}, onPartKey(DISFL(v-3).abbr);
            case 8, startAddError();
            case 9, onMultiKey();
        end
    end

    function syncPlayBtn()
        if isempty(btnBarPlay) || ~ishghandle(btnBarPlay), return; end
        if isPlaying, set(btnBarPlay,'String','Pause'); else, set(btnBarPlay,'String','Play'); end
    end

    function updateInfo()
        if isempty(infoTxt) || ~ishghandle(infoTxt), return; end
        if ~hasSel()
            s = {sprintf('Cursor: %s',fmtTime(cursorTime)),'No selection - drag on the timeline.'};
        else
            s = {sprintf('Selection: %s - %s',fmtTime(selStart),fmtTime(selEnd)), ...
                 sprintf('Length: %.3f s',selEnd-selStart)};
            if ~isempty(currentEvent) && currentEvent <= numel(events)
                s{2} = [s{2} '   |   ' kindLabel(currentEvent)];
            end
        end
        set(infoTxt,'String',s);
    end

    function updateAnnotDisplay()
        if isempty(annotDisp) || ~ishghandle(annotDisp), return; end
        if ~isempty(pendingNew)
            L = {'NEW TRANSCRIPT', sprintf('%s - %s',fmtTime(pendingNew.a),fmtTime(pendingNew.b)), '', ...
                'Type what was said, e.g.  t-t-table', 'It is split into pieces you can', ...
                'right-click and mark.', '', 'Enter = add,  Esc = cancel'};
        elseif isempty(currentEvent) || currentEvent > numel(events)
            L = {'Nothing selected.','', ...
                'Highlight speech -> double-click -> type it.', ...
                'Highlight a word -> Enter = fluent event.', ...
                'Right-click (two-finger click) a piece', 'to mark it fluent / disfluent / error.'};
        else
            k = currentEvent; ev = events(k); L = {};
            switch ev.kind
                case 'trans'
                    if any(strcmp({events(attachedTo(k)).kind},'fluent'))
                        L{end+1} = sprintf('FLUENT EVENT  "%s"',ev.transcript);
                    else
                        L{end+1} = sprintf('PIECE  "%s"',ev.transcript);
                    end
                case 'fluent', L{end+1} = sprintf('FLUENT WORD  "%s"',ev.transcript);
                case 'disfl',  L{end+1} = sprintf('DISFLUENCY:  %s  (%s)',ev.cat,disflName(ev.cat));
                otherwise,     L{end+1} = sprintf('ERROR:  %s  (%s)',ev.cat,errName(ev.cat));
            end
            L{end+1} = sprintf('%s - %s   (%.3f s)',fmtTime(ev.start),fmtTime(ev.end),ev.end-ev.start);
            [a,b] = currentTarget();
            L{end+1} = ['Marked as: ' pieceTypes(a,b)];
            L{end+1} = '';
            if strcmp(ev.kind,'disfl')
                L{end+1} = sprintf('What was said: "%s"',ev.transcript);
                f = findUid(ev.link);
                if isempty(f), L{end+1} = 'Resolves to: (no fluent word yet)';
                else, L{end+1} = sprintf('Resolves to: "%s"',events(f).transcript); end
                L{end+1} = '';
            end
            L{end+1} = 'Notes:  (N to edit)';
            L{end+1} = orNone(ev.notes,'(none)');
            u = groupUid(k);
            if u > 0, L = [groupLines(u) L]; end
        end
        if groupU > 0 && isempty(pendingNew)
            f = findUid(groupU);
            if ~isempty(f)
                L = [{sprintf('>> FLUENT EVENT SELECTED: "%s"',events(f).transcript), ...
                    '   Highlight a disfluent part + Enter', ...
                    '   W = whole episode,  Esc = done', ...
                    '   New speech: highlight + double-click',''} L];
            end
        end
        set(annotDisp,'String',L);
    end

    function L = groupLines(u)
        % e.g.   t- (rep)   t- (rep)   table (fluent)
        L = {};
        f = findUid(u); if isempty(f), return; end
        m = [f linkedAll(u)];
        [~,o] = sort([events(m).start]); m = m(o);
        L{end+1} = sprintf('"%s" BREAKDOWN',events(f).transcript);
        for q = m
            ev = events(q);
            switch ev.kind
                case 'fluent', lab = 'fluent'; tx = ev.transcript;
                case 'disfl',  lab = ev.cat;   tx = ev.transcript;
                otherwise,     lab = ev.cat;   tx = ['(' errName(ev.cat) ')'];
            end
            mark = '  '; if isequal(q,currentEvent), mark = '> '; end
            L{end+1} = sprintf('%s%s  %-6s %s',mark,fmtTime(ev.start),lab,tx); %#ok<AGROW>
        end
        L{end+1} = '';
    end

    function drawOverview()
        deleteValid(ovDyn); ovDyn = gobjects(0); ovCur = gobjects(0);
        if isempty(vid), return; end
        pv = patch(axOverview,'XData',[viewStart viewEnd viewEnd viewStart],'YData',[0 0 1 1], ...
            'FaceColor',OV_COLOR,'FaceAlpha',0.25,'EdgeColor',OV_COLOR,'LineWidth',1.5);
        ovDyn = [ovDyn pv];
        if ~isempty(events)
            kinds = {events.kind}; cats = {events.cat};
            for q = 1:numel(DISFL)
                ev = events(strcmp(kinds,'disfl') & strcmp(cats,DISFL(q).abbr));
                if isempty(ev), continue; end
                xs = reshape([[ev.start]; [ev.end]; nan(1,numel(ev))],1,[]);
                pe = plot(axOverview,xs,0.12*ones(size(xs)),'-','Color',DISFL(q).color,'LineWidth',4);
                ovDyn = [ovDyn pe]; %#ok<AGROW>
            end
            ev = events(strcmp(kinds,'error'));
            if ~isempty(ev)
                xs = reshape([[ev.start]; [ev.end]; nan(1,numel(ev))],1,[]);
                pe = plot(axOverview,xs,0.88*ones(size(xs)),'-','Color',ERR_COLOR,'LineWidth',4);
                ovDyn = [ovDyn pe];
            end
        end
        ovDyn = [ovDyn drawHandles(axOverview,viewStart,viewEnd)];
        ovCur = plot(axOverview,[cursorTime cursorTime],[0 1],'-','Color',CUR_COLOR,'LineWidth',1.2);
        ovDyn = [ovDyn ovCur];
        set(ovDyn,'PickableParts','none','HitTest','off');
    end

    function selectWholeView()
        selStart = viewStart; selEnd = viewEnd; currentEvent = [];
        cursorTime = selStart;
        redrawEvents(); updateSelectionGraphics(); updateCursorGraphics(); showFrameAt(cursorTime);
        setStatus(selHint(),'ok');
    end

    function s = selHint()
        if groupU > 0
            s = sprintf(['Highlighted %s - %s. Enter = disfluency / error type(s) | ' ...
                'double-click = new transcript | Esc = done with this word.'],fmtTime(selStart),fmtTime(selEnd));
        else
            s = sprintf(['Highlighted %s - %s (%.3f s). Double-click = type transcript | Enter = fluent ' ...
                'word (if on the transcript) | Space = play | right-click = menu.'], ...
                fmtTime(selStart),fmtTime(selEnd),selEnd-selStart);
        end
    end

    % ---------------- annotation menu ---------------------------------------

    % ---------------- mouse on the timelines --------------------------------
    % Press on a black edge handle -> resize.  Press inside an item -> drag to
    % move it (a transcript carries the annotations inside it along).  Press
    % elsewhere (or Shift + press) -> drag a new selection.
    function onAxDown(ax,panel)
        if nargin < 2, panel = ''; end
        if isempty(vid), return; end
        selType = get(fig,'SelectionType');
        cp = get(ax,'CurrentPoint'); t = clampT(cp(1,1)); y = cp(1,2);
        if strcmp(selType,'alt')
            % right-click: choose what the menu acts on, then MATLAB opens it
            commitEdit();
            rightClickPick(t,panel,y);
            return;
        end
        if strcmp(selType,'open')
            % double-click: transcript / notes for a highlight, or edit an item
            onDoubleClick(t,panel,y);
            return;
        end
        commitEdit();
        dragY = y;
        shift = any(strcmp(get(fig,'CurrentModifier'),'shift'));
        dragPanel = panel; resizeUndoPushed = false; movePushed = false; moveIdx = [];
        moveKids = []; moveKidsOrig = zeros(2,0); dragPrevCur = currentEvent; resizeNbr = [];
        dragMode = '';
        if ~shift, dragMode = hitEdge(ax,t); end
        if ~isempty(dragMode)
            if ~isempty(currentEvent)
                pushUndo(); resizeUndoPushed = true;
                if strcmp(events(currentEvent).kind,'trans')
                    resizeNbr = touchingPiece(currentEvent,dragMode);
                end
            end
        else
            hit = eventAtTime(t,panel,y);
            if (isempty(panel) || strcmp(panel,'transcript')) && ~isempty(hit)
                dragMode = 'select'; moveIdx = hit;   % drag = highlight, click = select
            elseif ~shift && ~isempty(hit)
                dragMode = 'move'; moveIdx = hit;
                moveOrig = [events(hit).start events(hit).end];
                if strcmp(events(hit).kind,'trans')
                    moveKids = attachedTo(hit);
                    moveKidsOrig = [[events(moveKids).start]; [events(moveKids).end]];
                end
            else
                dragMode = 'select';
            end
        end
        dragAxes = ax; dragStartT = t; downPix = get(fig,'CurrentPoint'); didDrag = false;
        set(fig,'WindowButtonMotionFcn',@(~,~)cb(@onDrag),'WindowButtonUpFcn',@(~,~)cb(@onUp));
    end

    function onDoubleClick(t,panel,y)
        % Double-click INSIDE the highlight (also when a rep / prol / td / error
        % on it is selected):
        %   spectrogram / waveform / Transcript strip -> type its transcript
        %     (or edit the transcript piece already there)
        %   Notes strip -> notes for the selected item, else the transcript
        %     piece there, else the item under the pointer
        % Double-click an item outside the highlight: edit its text / notes.
        commitEdit();
        inHighlight = hasSel() && t >= selStart && t <= selEnd;
        if inHighlight
            if strcmp(panel,'notes')
                if ~isempty(currentEvent) && currentEvent <= numel(events)
                    k = currentEvent;
                else
                    k = pieceUnder(selStart,selEnd,t);
                    if isempty(k), k = eventAtTime(t,'notes',y); end
                end
                if isempty(k)
                    setStatus(['Nothing here to attach notes to yet - double-click the highlight ' ...
                        'on the spectrogram or Transcript strip to type its transcript first.'],'warn');
                    return;
                end
                startEditing(k,'notes'); return;
            end
            if groupU > 0, exitGroup(); end
            if isempty(overlapIdx(selStart,selEnd,{'trans'},[]))
                beginNewText('trans');                                % new transcript
            else
                startEditing(pieceUnder(selStart,selEnd,t),'text');   % edit existing piece
            end
            return;
        end
        idx = eventAtTime(t,panel,y);
        if isempty(idx) || isempty(panel), return; end
        if strcmp(panel,'notes'), startEditing(idx,'notes'); else, startEditing(idx,'text'); end
    end

    function k = pieceUnder(a,b,t)
        % Transcript piece under time t, else the one overlapping [a b] most.
        k = [];
        if isempty(events), return; end
        tr = find(strcmp({events.kind},'trans'));
        if isempty(tr), return; end
        at = tr([events(tr).start] <= t & [events(tr).end] >= t);
        if ~isempty(at), k = at(1); return; end
        ov = min(b,[events(tr).end]) - max(a,[events(tr).start]);
        [m,i] = max(ov);
        if m > TOL, k = tr(i); end
    end

    function onDrag()
        cp = get(dragAxes,'CurrentPoint'); t = clampT(cp(1,1));
        if norm(get(fig,'CurrentPoint')-downPix) > 3, didDrag = true; end
        if ~didDrag, return; end
        switch dragMode
            case 'left'
                lo = 0; if ~isempty(resizeNbr), lo = events(resizeNbr).start + MINSEL; end
                selStart = max(lo, min(t, selEnd-MINSEL));
                if ~isempty(currentEvent)
                    events(currentEvent).start = selStart; clearGlobal(currentEvent);
                    if strcmp(events(currentEvent).kind,'trans'), syncTok(currentEvent);
                    else, events(currentEvent).tok = 0; end
                    if ~isempty(resizeNbr)
                        events(resizeNbr).end = selStart; clearGlobal(resizeNbr); syncTok(resizeNbr);
                    end
                    redrawEvents();
                end
            case 'right'
                hi = dur; if ~isempty(resizeNbr), hi = events(resizeNbr).end - MINSEL; end
                selEnd = min(hi, max(t, selStart+MINSEL));
                if ~isempty(currentEvent)
                    events(currentEvent).end = selEnd; clearGlobal(currentEvent);
                    if strcmp(events(currentEvent).kind,'trans'), syncTok(currentEvent);
                    else, events(currentEvent).tok = 0; end
                    if ~isempty(resizeNbr)
                        events(resizeNbr).start = selEnd; clearGlobal(resizeNbr); syncTok(resizeNbr);
                    end
                    redrawEvents();
                end
            case 'move'
                if ~isempty(moveIdx) && ~movePushed
                    pushUndo(); movePushed = true;
                    currentEvent = moveIdx;
                end
                a = moveOrig(1) + (t - dragStartT); b = moveOrig(2) + (t - dragStartT);
                if a < 0,   b = b - a; a = 0; end
                if b > dur, a = a - (b-dur); b = dur; end
                d = a - moveOrig(1);
                selStart = a; selEnd = b;
                if ~isempty(moveIdx)
                    events(moveIdx).start = a; events(moveIdx).end = b; clearGlobal(moveIdx);
                    if ~strcmp(events(moveIdx).kind,'trans'), events(moveIdx).tok = 0; end
                    for j = 1:numel(moveKids)
                        kk = moveKids(j);
                        events(kk).start = moveKidsOrig(1,j) + d;
                        events(kk).end   = moveKidsOrig(2,j) + d;
                        clearGlobal(kk);
                    end
                    redrawEvents();
                end
            otherwise
                selStart = min(dragStartT,t); selEnd = max(dragStartT,t); currentEvent = []; moveIdx = [];
        end
        updateSelectionGraphics();
        if ~isnan(selStart), syncTimeFields(selStart,selEnd); end
        drawnow limitrate;
    end

    function onUp()
        set(fig,'WindowButtonMotionFcn',@(~,~)onHover(),'WindowButtonUpFcn','');
        m0 = dragMode; dragMode = ''; panel = dragPanel; t = dragStartT;
        if any(strcmp(m0,{'left','right'}))
            if didDrag
                cursorTime = selStart;
                if ~isempty(currentEvent)
                    if ~partitionOK(currentEvent)
                        revertLast('Not allowed: overlap - resize undone.');
                    else
                        relinkAll();
                        dirty = true; updateTitle(); updateEventList(); redrawEvents();
                        setStatus(sprintf('Resized %s.',eventLabel(currentEvent)),'ok');
                    end
                else
                    setStatus(selHint(),'info');
                end
            elseif resizeUndoPushed && ~isempty(undoStack)
                undoStack(end) = [];
            end
        elseif strcmp(m0,'move') && didDrag
            cursorTime = selStart;
            if ~isempty(moveIdx)
                if ~partitionOK(moveIdx)
                    revertLast('Not allowed: overlap - move undone.');
                else
                    relinkAll();
                    dirty = true; updateTitle(); updateEventList(); redrawEvents();
                    setStatus(sprintf('Moved %s (Ctrl+Z to undo).',eventLabel(moveIdx)),'ok');
                end
            end
        elseif didDrag
            cursorTime = selStart; redrawEvents();
            setStatus(selHint(),'info');
        else
            % plain click: select what is under the pointer
            idx = moveIdx; if isempty(idx), idx = eventAtTime(t,panel,dragY); end
            if ~isempty(idx)
                currentEvent = idx; selStart = events(idx).start; selEnd = events(idx).end;
                cursorTime = selStart; syncListSelection();
                setStatus(sprintf(['Selected %s. Right-click (two-finger click) for options, ' ...
                    'double-click to edit its text, Space to play.'],eventLabel(idx)),'info');
            elseif (isempty(panel) || strcmp(panel,'transcript')) && hasSel() && isempty(currentEvent) ...
                    && t > selStart && t < selEnd
                cursorTime = t;                       % click inside the highlight: play point
                setStatus(sprintf(['Play point %s - Space plays from here. Double-click the ' ...
                    'highlight to type its transcript.'],fmtTime(t)),'info');
            elseif isempty(panel)
                selStart = NaN; selEnd = NaN; currentEvent = []; cursorTime = t;
            end
            redrawEvents();
        end
        moveIdx = []; moveKids = []; moveKidsOrig = zeros(2,0); resizeNbr = [];
        updateSelectionGraphics(); updateCursorGraphics(); showFrameAt(cursorTime);
    end

    function onHover()
        try
            if isempty(vid) || isEditing || isempty(fig) || ~ishghandle(fig), return; end
            ptr = 'arrow';
            axs = [axSpec axWave axTrans axDE axNotes];
            pnl = {'','','transcript','disferr','notes'};
            for q = 1:numel(axs)
                hx = axs(q);
                cp = get(hx,'CurrentPoint'); xl = get(hx,'XLim'); yl = get(hx,'YLim');
                if cp(1,1)>=xl(1) && cp(1,1)<=xl(2) && cp(1,2)>=yl(1) && cp(1,2)<=yl(2)
                    x = cp(1,1);
                    if ~isempty(hitEdge(hx,x))
                        ptr = 'left';
                    elseif ~isempty(eventAtTime(x,pnl{q},cp(1,2)))
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
            updateFrameStride();
            updateSelectionGraphics(); positionBoxes();
        catch
        end
    end

    % ---------------- keyboard ----------------------------------------------
    function onKey(e)
        k = e.Key;
        if strcmp(k,'f1'), showHelp(); return; end
        if isempty(vid)
            if any(strcmp(k,{'space','r','b','p','t','e','f','m','w','n','return'}))
                setStatus('Load a video first (click "Load video").','warn');
            end
            return;
        end
        if isEditing
            if strcmp(k,'escape'), cancelTyping(); end
            % backup in case the box's own key handler did not see Enter; the
            % box callback (which has the up-to-date text) does the saving
            if any(strcmp(k,{'return','enter'})), lastBoxKey = 'return'; end
            return;
        end
        ctrl  = any(strcmp(e.Modifier,'control')) || any(strcmp(e.Modifier,'command'));
        shift = any(strcmp(e.Modifier,'shift'));
        if ctrl
            switch k
                case 'o', zoomAbout(cursorCenter(),2);
                case 'i', zoomAbout(cursorCenter(),0.5);
                case 'n', zoomToSelection();
                case 'a', fullView();
                case 's', saveAnnotations();
                case 'z', undo();
            end
            return;
        end
        switch k
            case 'space',                togglePlay();
            case 'escape'
                if justCommitted(), return; end      % Esc already used by a text box
                if isPlaying, stopPlay(); setStatus('Stopped.','info');
                elseif groupU > 0, exitGroup();
                else, clearSelection(); end
            case {'return','enter'}
                if justCommitted(), return; end      % Enter already used by a text box
                onEnter();
            case {'delete','backspace'}, deleteCurrentEvent();
            case 'leftarrow',            stepCursor(-1,shift);
            case 'rightarrow',           stepCursor(1,shift);
            case ERR_KEY,                startAddError();
            case FLU_KEY,                onPartKey('flu');
            case 'm',                    onMultiKey();
            case 'w',                    wholeEpisode();
            case 'n',                    editNotes();
            case 'l',                    setLoop(~loopOn);
            case {'1','2','3','4'},      setSpeedIdx(str2double(k));
            otherwise
                q = find(strcmp({DISFL.key},k),1);
                if ~isempty(q), onPartKey(DISFL(q).abbr); end
        end
    end

    function boxKey(e)
        lastBoxKey = e.Key;
        if strcmp(e.Key,'escape'), cancelTyping(); end
    end

    function insertSymbol(src)
        % IPA button: type the symbol into the box being edited (else copy it).
        sym = get(src,'String'); if iscell(sym), sym = sym{1}; end
        box = typingBox;
        if isEditing && ~isempty(box) && ishghandle(box) && ...
                ~any(strcmp(get(box,'Tag'),{'tfrom','tto'}))
            str = get(box,'String'); if iscell(str), str = strjoin(str,' '); end
            set(box,'String',[char(str) sym],'Enable','on');
            uicontrol(box);
            setStatus(sprintf('Inserted %s - keep typing; Enter saves.',sym),'info');
        else
            dropFocus(src);
            try, clipboard('copy',sym); catch, end
            setStatus(sprintf(['%s copied. Start typing in a transcript / notes box first ' ...
                'to insert it directly.'],sym),'info');
        end
    end

    % ---------------- adding items ------------------------------------------
    % ---------------- Enter workflow: fluent word -> its disfluent parts ------
    % 0. Highlight new speech + double-click -> type the transcript.
    % 1. Highlight a word in the transcript (or click its piece) + Enter
    %    -> "fluent event selected" (group mode).
    % 2. Highlight each disfluent part (t-, t-t, a syllable, or W = the whole
    %    episode) + Enter -> tick one or more disfluency / error types
    %    (keys R P B T, 1-6) -> Enter. Repeat. Esc = done.
    % Enter on a highlight with no transcript under it (outside group mode)
    % = type the transcript too.
    function onEnter()
        if groupU > 0 && isempty(findUid(groupU)), groupU = 0; end
        [a,b,tk] = currentTarget();
        if isnan(a)
            if groupU > 0
                setStatus('Highlight a disfluent part (e.g. t-) and press Enter - or Esc when done.','info');
            else
                setStatus('Drag to highlight first: a word in the transcript, or new speech to transcribe.','warn');
            end
            return;
        end
        if groupU > 0
            f = findUid(groupU);
            if abs(a-events(f).start) < 1e-3 && abs(b-events(f).end) < 1e-3
                setStatus(['That is the fluent word itself. Highlight a disfluent part (e.g. t-) and ' ...
                    'press Enter, W for the whole episode, or Esc when done.'],'info');
                return;
            end
            openEventPicker(); return;
        end
        % not in group mode
        if isempty(overlapIdx(a,b,{'trans'},[]))
            % no transcript here yet - even if you already tagged the whole
            % episode R / P / B / T: Enter types the transcript
            if ~hasSel(), selStart = a; selEnd = b; end
            beginNewText('trans'); return;
        end
        if ~isempty(currentEvent) && currentEvent <= numel(events) && ...
                any(strcmp(events(currentEvent).kind,{'disfl','error'}))
            openEventPicker(); return;                 % change the types of this part
        end
        fluentFromTarget();
    end

    function beginNewText(kind)
        % kind 'trans': type the transcript for the selection
        if ~hasSel()
            msgbox('Drag across the spectrogram / waveform to select first.','Select first','help','replace');
            return;
        end
        a = selStart; b = selEnd;
        ov = overlapIdx(a,b,{'trans'},[]);
        if ~isempty(ov)
            setStatus(sprintf(['That overlaps the transcript "%s" - click it on the Transcript ' ...
                'strip to edit it, or trim the selection.'],events(ov(1)).transcript),'warn');
            beep; return;
        end
        if isEditing, commitEdit(); end
        groupU = 0;                                   % a new transcript ends the current word
        currentEvent = [];
        pendingNew = struct('kind',kind,'a',a,'b',b);
        box = boxFor(kind);
        set(box,'String','','UserData',-1);
        typingBox = box; isEditing = true; lastBoxKey = '';
        redrawEvents();                               % also positions the box
        set(box,'Enable','on'); uicontrol(box);
        setStatus('TRANSCRIPT: type what was said, Enter to add it, Esc to cancel.','info');
    end

    function onPartKey(code)
        % F = fluent word (same as Enter); R / P / B / T = add that disfluency
        % straight away (in group mode it joins the selected fluent word)
        if strcmp(code,'flu'), fluentFromTarget(); return; end
        [a,b,tk] = currentTarget();
        if isnan(a), setStatus('Highlight the part first.','warn'); return; end
        [a,b,tk] = snapTarget(a,b,tk);
        selStart = a; selEnd = b;
        if ~isempty(tk), currentEvent = tk; elseif ~isempty(currentEvent) && ...
                ~any(strcmp(events(currentEvent).kind,{'disfl','error'})), currentEvent = []; end
        addAnnotation('disfl',code);
    end

    function w = normWordKeep(tok)
        % the word as shown, minus surrounding punctuation (case kept)
        w = regexprep(strtrim(char(tok)),'^[\.,;:!?"''()\[\]-]+|[\.,;:!?"''()\[\]-]+$','');
        if isempty(w), w = char(tok); end
    end

    function wi = wordPicker(sk,a,b)
        % Ask which word of transcript sk the fluent piece [a b] is.
        wi = [];
        toks = regexp(events(sk).transcript,'\S+','match');
        if isempty(toks)
            msgbox('That transcript is empty - type it first (click it on the Transcript strip).', ...
                'Fluent word','help','replace');
            return;
        end
        if numel(toks) == 1, wi = 1; return; end
        n = numel(toks);
        frac = ((a+b)/2 - events(sk).start)/max(events(sk).end - events(sk).start,eps);
        est = min(max(ceil(frac*n),1),n);
        used = [events(fluentsInSeg(sk)).wi];
        cand = setdiff(1:n,used); if isempty(cand), cand = 1:n; end
        [~,m] = min(abs(cand-est)); def = cand(m);
        names = arrayfun(@(i)sprintf('%d.  %s',i,toks{i}),1:n,'UniformOutput',false);
        [sel,ok] = listdlg('PromptString',{'Which word is this fluent piece?','(Enter = OK, Esc = cancel)'}, ...
            'SelectionMode','single','ListString',names,'InitialValue',def, ...
            'Name','Fluent word','ListSize',[240 min(320,22*n+30)]);
        if ok, wi = sel; end
    end

    function promptPartText(k)
        % Optional: type what was actually said ('t', 's'...). Blocks: no speech.
        f = findUid(events(k).link);
        if isempty(f), tail = ' (not connected to a fluent word yet)';
        else, tail = sprintf(' -> "%s"',events(f).transcript); end
        if DISFL_HAS_TRANSCRIPT && ~strcmp(events(k).cat,'bl')
            startEditing(k,'text');
            setStatus(sprintf(['Added %s%s. Type what was said (e.g. "t") if you like - Enter saves, ' ...
                'Esc skips. Click another type to add a second event.'],events(k).cat,tail),'info');
        else
            setStatus(sprintf('Added %s%s. Ctrl+Z to undo.',events(k).cat,tail),'ok');
        end
    end

    function setPartType(k,code)
        % Change type button: retype a fluent word / disfluency.
        ev = events(k);
        if (strcmp(code,'flu') && strcmp(ev.kind,'fluent')) || strcmp(ev.cat,code)
            setStatus(sprintf('%s is already that type.',kindLabel(k)),'info'); return;
        end
        if strcmp(code,'flu')
            if ~isempty(overlapIdx(ev.start,ev.end,{'fluent'},k))
                setStatus('That overlaps another fluent word.','warn'); return;
            end
            sk = segOf(k); wi = 0; word = ev.transcript;
            if ~isempty(sk)
                wi = wordPicker(sk,ev.start,ev.end); if isempty(wi), return; end
                toks = regexp(events(sk).transcript,'\S+','match'); word = normWordKeep(toks{wi});
            end
            pushUndo();
            events(k).kind = 'fluent'; events(k).cat = ''; events(k).link = 0;
            events(k).wi = wi; events(k).transcript = word;
        else
            pushUndo();
            if strcmp(ev.kind,'fluent'), events(k).transcript = ''; events(k).wi = 0; end
            events(k).kind = 'disfl'; events(k).cat = code;
        end
        relinkAll();
        eventsChanged();
        setStatus(sprintf('Changed to %s (Ctrl+Z to undo).',kindLabel(k)),'ok');
    end

    function startAddError()
        [a,~] = currentTarget();
        if isnan(a)
            setStatus('Highlight some speech or click a transcript piece first.','warn'); return;
        end
        q = errorPicker('Add error');
        if isempty(q), setStatus('Cancelled - no error added.','info'); return; end
        addAnnotation('error',ERRS(q).abbr);
    end

    function k = newEvent(a,b,kind,cat,tr,link)
        s = struct('start',a,'end',b,'kind',kind,'cat',cat,'transcript',tr,'notes','', ...
            'uid',nextUid,'link',link,'wi',0,'tok',0,'gStart',NaN,'gStop',NaN);
        nextUid = nextUid + 1;
        events(end+1) = s; k = numel(events);
    end

    function idx = eventAtTime(t,panel,y)
        % Shortest item at time t among the kinds shown where you clicked.
        idx = [];
        if isempty(events), return; end
        if nargin < 2, panel = ''; end
        if nargin < 3, y = NaN; end
        hit = find([events.start] <= t & [events.end] >= t);
        if isempty(hit), return; end
        kinds = {events(hit).kind};
        switch panel
            case 'transcript'
                pref = hit(strcmp(kinds,'trans'));
            case 'disferr'
                cand = hit(~ismember(kinds,{'trans','fluent'})); pref = [];
                if ~isnan(y)
                    for c = cand
                        bnd = laneBand(c);
                        if y >= bnd(1)-0.01 && y <= bnd(2)+0.01, pref(end+1) = c; end %#ok<AGROW>
                    end
                end
                if isempty(pref), pref = cand; end
            case 'notes'
                pref = hit;
            otherwise
                pref = hit(~ismember(kinds,{'trans','fluent'}));
                if isempty(pref), pref = hit(strcmp(kinds,'trans')); end   % the word piece
        end
        if isempty(pref), return; end
        [~,m] = min([events(pref).end] - [events(pref).start]);
        idx = pref(m);
    end

    function idx = overlapIdx(a,b,kinds,exclude)
        idx = [];
        if isempty(events), return; end
        idx = find([events.start] < b - TOL & [events.end] > a + TOL & ...
            ismember({events.kind},kinds));
        idx = setdiff(idx,exclude);
    end

    function sk = segContaining(a,b)
        sk = [];
        if isempty(events), return; end
        sk = find(strcmp({events.kind},'trans') & [events.start] <= a + TOL & ...
            [events.end] >= b - TOL,1);
    end

    function sk = segOf(k)
        sk = [];
        if isempty(k) || k > numel(events) || strcmp(events(k).kind,'trans'), return; end
        sk = segContaining(events(k).start,events(k).end);
    end

    function idx = insideSeg(sk)
        idx = find(~strcmp({events.kind},'trans') & [events.start] >= events(sk).start - TOL & ...
            [events.end] <= events(sk).end + TOL);
    end

    function idx = fluentsInSeg(sk)
        idx = [];
        if isempty(sk), return; end
        idx = find(strcmp({events.kind},'fluent') & [events.start] >= events(sk).start - TOL & ...
            [events.end] <= events(sk).end + TOL);
        [~,o] = sort([events(idx).start]); idx = idx(o);
    end

    function nf = nextFluent(sk,t)
        % First fluent word in transcript sk that starts at or after time t.
        nf = [];
        fl = fluentsInSeg(sk);
        fl = fl([events(fl).start] >= t - TOL);
        if ~isempty(fl), nf = fl(1); end
    end

    function sp = groupSpan(u)
        m = [findUid(u) linkedAll(u)];
        sp = [min([events(m).start]) max([events(m).end])];
    end

    function ok = partitionOK(k)
        % Only rules: transcript pieces don't overlap each other, and fluent
        % words don't overlap each other. Everything else may overlap.
        ok = true;
        if isempty(k) || k > numel(events), return; end
        ev = events(k);
        switch ev.kind
            case 'trans',  ok = isempty(overlapIdx(ev.start,ev.end,{'trans'},k));
            case 'fluent', ok = isempty(overlapIdx(ev.start,ev.end,{'fluent'},k));
        end
    end

    function idx = linkedAll(u)
        % Disfluencies and errors connected to fluent word uid u.
        idx = [];
        if isempty(events) || u <= 0, return; end
        idx = find(ismember({events.kind},{'disfl','error'}) & [events.link] == u);
    end

    function u = groupUid(k)
        % uid of the fluent word item k belongs to (0 = none). For a transcript
        % piece: the group of the events on it.
        u = 0;
        if isempty(k) || k > numel(events), return; end
        switch events(k).kind
            case 'fluent', u = events(k).uid;
            case {'disfl','error'}
                if events(k).link > 0 && ~isempty(findUid(events(k).link)), u = events(k).link; end
            case 'trans'
                for j = attachedTo(k)
                    u = groupUid(j); if u > 0, return; end
                end
        end
    end

    % ---------------- mouse on the timelines --------------------------------
    % Left-drag on the spectrogram / waveform = highlight. Click an item to
    % select it; drag it on a strip to move it; drag black edges to resize
    % (between two transcript pieces both edges move together). Right-click
    % (two-finger click) = annotation menu. Double-click = edit text.

    function n = fluentNumber(k)
        fl = find(strcmp({events.kind},'fluent'));
        [~,o] = sort([events(fl).start]); fl = fl(o);
        n = find(fl == k,1); if isempty(n), n = 0; end
    end

    function deleteCurrentEvent()
        if isempty(currentEvent) || currentEvent > numel(events)
            setStatus('Nothing selected - click an item first, then Delete.','warn'); return;
        end
        commitEditQuiet();
        k = currentEvent; ev = events(k);
        if strcmp(ev.kind,'trans')
            att = attachedTo(k);
            if ~isempty(att) && all(strcmp({events(att).kind},'fluent'))
                pushUndo(); events([k att]) = [];
            elseif ~isempty(att)
                c = questdlg(sprintf('Delete the piece "%s" and what it is marked as (%s)?', ...
                    ev.transcript,pieceTypes(ev.start,ev.end)),'Delete piece', ...
                    'Delete all','Piece only','Cancel','Delete all');
                switch c
                    case 'Delete all', pushUndo(); events([k att]) = [];
                    case 'Piece only', pushUndo(); [events(att).tok] = deal(0); events(k) = [];
                    otherwise, return;
                end
            else
                pushUndo(); events(k) = [];
            end
        else
            pushUndo(); events(k) = [];
        end
        relinkAll();
        if groupU > 0 && isempty(findUid(groupU)), groupU = 0; end
        currentEvent = []; selStart = NaN; selEnd = NaN;
        eventsChanged(); updateSelectionGraphics();
        setStatus('Deleted (Ctrl+Z to undo).','ok');
    end

    function changeTypeDialog()
        if isempty(currentEvent)
            msgbox('Select an item first (click it on the timeline or in the Events list).', ...
                'Change type','help','replace'); return;
        end
        k = currentEvent;
        switch events(k).kind
            case 'trans'
                msgbox('A transcript has no type. Select a fluent word, disfluency or error.', ...
                    'Change type','help','replace');
            case {'fluent','disfl'}
                codes = [{'flu'} {DISFL.abbr}];
                names = [{'flu - fluent word'} arrayfun(@(d)sprintf('%s - %s',d.abbr,d.name),DISFL,'UniformOutput',false)];
                cur = find(strcmp(codes,events(k).cat),1);
                if strcmp(events(k).kind,'fluent') || isempty(cur), cur = 1; end
                [sel,ok] = listdlg('PromptString','Choose the type:','SelectionMode','single', ...
                    'ListString',names,'InitialValue',cur,'Name','Change type','ListSize',[260 110]);
                if ok, setPartType(k,codes{sel}); end
            otherwise
                q = errorPicker('Change error type');
                if ~isempty(q) && ~strcmp(events(k).cat,ERRS(q).abbr)
                    pushUndo(); events(k).cat = ERRS(q).abbr; eventsChanged();
                    setStatus(sprintf('Changed to %s (Ctrl+Z to undo).',ERRS(q).abbr),'ok');
                end
        end
    end

    function editCurrentEvent()
        if isempty(currentEvent)
            msgbox('Select an item first (click it on the timeline or in the Events list).', ...
                'Edit text','help','replace'); return;
        end
        k = currentEvent;
        if events(k).start < viewStart || events(k).end > viewEnd, selectEvent(k,true); end
        if strcmp(events(k).kind,'error'), startEditing(k,'notes');
        else, startEditing(k,'text'); end
    end

    % ---------------- events list -------------------------------------------
    function updateEventList()
        listMap = []; listKind = {};
        if isempty(events)
            set(lstEvents,'Data',cell(0,5),'BackgroundColor',[1 1 1]); return;
        end
        fl = find(strcmp({events.kind},'fluent'));
        [~,o] = sort([events(fl).start]); fl = fl(o);
        expandedU = expandedU(ismember(expandedU,[events(fl).uid]));
        selU = 0;
        if ~isempty(currentEvent) && currentEvent <= numel(events), selU = groupUid(currentEvent); end
        if groupU > 0, selU = groupU; end
        data = cell(0,5); bgc = zeros(0,3);
        for f = fl
            u = events(f).uid; kids = linkedAll(u);
            [~,o] = sort([events(kids).start]); kids = kids(o);
            isOpen = any(expandedU == u);
            if isempty(kids), arrow = ' ';
            elseif isOpen, arrow = char(9662);              % down triangle
            else, arrow = char(9656); end                    % right triangle
            lab = events(f).transcript;
            if ~isempty(kids), lab = sprintf('%s  (%d)',lab,numel(kids)); end
            data(end+1,:) = {arrow, fmtTime(events(f).start), fmtTime(events(f).end), ...
                lab, oneLine(fluentNotes(f))}; %#ok<AGROW>
            if u == selU, bgc(end+1,:) = [0.80 0.88 1.00]; else, bgc(end+1,:) = [1 1 1]; end %#ok<AGROW>
            listMap(end+1) = f; listKind{end+1} = 'fluent'; %#ok<AGROW>
            if isOpen
                for j = kids
                    data(end+1,:) = {'', fmtTime(events(j).start), fmtTime(events(j).end), ...
                        ['   ' eventCode(j) '  ' events(j).cat], oneLine(events(j).notes)}; %#ok<AGROW>
                    bgc(end+1,:) = 0.55 + 0.45*eventColor(j); %#ok<AGROW>
                    listMap(end+1) = j; listKind{end+1} = 'child'; %#ok<AGROW>
                end
            end
        end
        % disfluencies / errors not connected to a fluent word
        orph = find(ismember({events.kind},{'disfl','error'}) & ...
            arrayfun(@(e) e.link <= 0 || isempty(findUid(e.link)),events));
        if ~isempty(orph)
            [~,o] = sort([events(orph).start]); orph = orph(o);
            if unresolvedOpen, arrow = char(9662); else, arrow = char(9656); end
            data(end+1,:) = {arrow, '', '', sprintf('(unresolved)  (%d)',numel(orph)), ''};
            bgc(end+1,:) = [0.92 0.92 0.92]; listMap(end+1) = 0; listKind{end+1} = 'unresolved';
            if ~unresolvedOpen, orph = []; end
            for j = orph
                data(end+1,:) = {'', fmtTime(events(j).start), fmtTime(events(j).end), ...
                    ['   ' eventCode(j) '  ' events(j).cat], oneLine(events(j).notes)}; %#ok<AGROW>
                bgc(end+1,:) = 0.55 + 0.45*eventColor(j); %#ok<AGROW>
                listMap(end+1) = j; listKind{end+1} = 'child'; %#ok<AGROW>
            end
        end
        if isempty(data)
            data = {'', '', '', '(no fluent events yet)', ''}; bgc = [1 1 1];
            listMap = 0; listKind = {'header'};
        end
        set(lstEvents,'Data',data,'BackgroundColor',bgc);
    end

    function s = kindLabel(k)
        ev = events(k);
        switch ev.kind
            case 'trans',  s = sprintf('transcript "%s"',ev.transcript);
            case 'fluent', s = sprintf('fluent word "%s"',ev.transcript);
            case 'disfl',  s = sprintf('%s "%s"',ev.cat,ev.transcript);
            otherwise,     s = sprintf('%s error',ev.cat);
        end
    end

    % ---------------- text editing ------------------------------------------
    % A white box opens over the item. Enter saves, Esc cancels. Clicking an
    % IPA button while typing inserts the symbol and keeps the box open.
    function startEditing(idx,panel)
        if isempty(vid) || isempty(idx) || idx > numel(events), return; end
        if ~isequal(currentEvent,idx)
            currentEvent = idx; selStart = events(idx).start; selEnd = events(idx).end;
            syncListSelection(); updateSelectionGraphics();
        end
        redrawEvents();
        beginTyping(panel);
    end

    function box = boxFor(kind)
        if strcmp(kind,'trans'), box = transBox; else, box = deBox; end
    end

    function beginTyping(panel)
        if isempty(vid) || isempty(currentEvent) || currentEvent > numel(events), return; end
        if isEditing && ~isempty(typingBox) && ishghandle(typingBox)
            commitBox(typingBox);
        end
        k = currentEvent;
        if strcmp(panel,'notes') || strcmp(events(k).kind,'error')
            box = notesBox; fld = 'notes';
        else
            box = boxFor(events(k).kind); fld = 'transcript';
        end
        set(box,'String',events(k).(fld),'UserData',k);
        typingBox = box; isEditing = true; lastBoxKey = '';
        positionBoxes();
        if ~strcmp(get(box,'Visible'),'on'), typingBox = []; isEditing = false; return; end
        set(box,'Enable','on');
        uicontrol(box);
        if strcmp(fld,'notes'), what = 'notes'; else, what = 'text'; end
        setStatus(sprintf('Typing the %s for %s - Enter saves, Esc cancels.',what,eventLabel(k)),'info');
    end

    function commitBox(box,force)
        % force = false for the box's own callback: only Enter saves, so focus
        % moving to an IPA button keeps the box open.
        if nargin < 2, force = true; end
        if isempty(box) || ~ishghandle(box), return; end
        tag = get(box,'Tag');
        if any(strcmp(tag,{'tfrom','tto'}))
            if ~isequal(typingBox,box), set(box,'Enable','inactive'); return; end
            set(box,'Enable','inactive'); typingBox = []; isEditing = false;
            lastKeyCommit = now*86400;
            if isempty(strtrim(char(get(box,'String'))))
                set(box,'String',get(box,'UserData'));
            end
            applyTimeFields(tag); return;
        end
        if ~isequal(typingBox,box)          % stale callback after Esc / commit
            set(box,'Enable','inactive','Visible','off'); return;
        end
        if ~force && ~any(strcmp(lastBoxKey,{'return','enter'})), return; end
        lastBoxKey = '';
        idx = get(box,'UserData');
        if isequal(box,notesBox), field = 'notes'; else, field = 'transcript'; end
        str = get(box,'String');
        if iscell(str), str = strjoin(str,' ');
        elseif size(str,1) > 1, str = strjoin(cellstr(str),' '); end
        str = strtrim(str);
        set(box,'Enable','inactive','Visible','off');
        typingBox = []; isEditing = false; lastKeyCommit = now*86400;

        if isequal(idx,-1)                   % finishing a new transcript
            pn = pendingNew; pendingNew = [];
            if isempty(pn), redrawEvents(); updateSelBar(); return; end
            if isempty(str)
                redrawEvents(); updateSelBar();
                setStatus('Nothing added - it needs text. Double-click the highlight to try again.','warn');
                return;
            end
            pushUndo();
            addPieces(pn.a,pn.b,str);
            currentEvent = []; selStart = NaN; selEnd = NaN;
            eventsChanged(); updateSelectionGraphics();
            setStatus(sprintf(['Added "%s". Drag the edges between pieces to fit the audio, then ' ...
                'highlight the fluent word (e.g. table) and press Enter.'],str),'ok');
            return;
        end

        if isempty(idx) || idx > numel(events), redrawEvents(); return; end
        old = events(idx).(field);
        if ~strcmp(old,str)
            pushUndo();
            events(idx).(field) = str;
            if strcmp(field,'transcript') && strcmp(events(idx).kind,'trans')
                % events on this piece take its new text
                for j = attachedTo(idx)
                    if ~strcmp(events(j).kind,'error'), events(j).transcript = tokenWord(str); end
                end
            end
            dirty = true; updateTitle(); updateEventList(); redrawEvents();
            setStatus(sprintf('Saved %s for %s.',field,eventLabel(idx)),'ok');
        else
            redrawEvents();
        end
        if ~isempty(batchUids) && strcmp(field,'transcript') && ~isempty(str)
            % several events added at once (M): they share what was said
            for u = batchUids
                j = findUid(u);
                if ~isempty(j) && strcmp(events(j).kind,'disfl'), events(j).transcript = str; end
            end
            dirty = true; updateEventList(); redrawEvents();
        end
        batchUids = [];
        updateSelBar();
    end

    function remapWords(sk)
        % After a transcript is edited, keep each fluent word pointing at the
        % same word (nearest match), or unhighlight it if the word is gone.
        for f = fluentsInSeg(sk), remapOne(f); end
    end

    function remapOne(f)
        sk = segOf(f); if isempty(sk), return; end
        toks = cellfun(@normWord,regexp(events(sk).transcript,'\S+','match'),'UniformOutput',false);
        w = normWord(events(f).transcript); wi = events(f).wi;
        if wi >= 1 && wi <= numel(toks) && strcmp(toks{wi},w), return; end
        m = find(strcmp(toks,w));
        if isempty(m), events(f).wi = 0;
        else, [~,j] = min(abs(m - max(wi,1))); events(f).wi = m(j); end
    end

    function cancelTyping()
        box = typingBox;
        if isempty(box) || ~ishghandle(box), isEditing = false; typingBox = []; return; end
        typingBox = []; isEditing = false; lastKeyCommit = now*86400; lastBoxKey = '';
        batchUids = [];
        tag = get(box,'Tag');
        if any(strcmp(tag,{'tfrom','tto'}))
            set(box,'Enable','inactive','String',get(box,'UserData'));
            setStatus('Time entry cancelled.','info'); return;
        end
        set(box,'Enable','inactive','Visible','off');
        if isequal(get(box,'UserData'),-1)
            pendingNew = [];
            setStatus('Cancelled - nothing added.','info');
        else
            setStatus('Edit cancelled - text unchanged.','info');
        end
        redrawEvents(); updateSelBar();
    end

    function commitEditQuiet()
        isEditing = false; typingBox = []; pendingNew = [];
        for bx = [transBox deBox notesBox]
            if ~isempty(bx) && ishghandle(bx), set(bx,'Enable','inactive','Visible','off'); end
        end
        for bx = [edtFrom edtTo]
            if ~isempty(bx) && ishghandle(bx), set(bx,'Enable','inactive'); end
        end
    end

    function positionBoxes()
        % A white text box is shown only while you type in it.
        boxes = {transBox, axTrans; deBox, axDE; notesBox, axNotes};
        for q = 1:size(boxes,1)
            box = boxes{q,1};
            if isempty(box) || ~ishghandle(box), continue; end
            if ~isequal(typingBox,box), set(box,'Visible','off'); continue; end
            k = get(box,'UserData');
            if isequal(k,-1)
                if isempty(pendingNew), set(box,'Visible','off'); continue; end
                a = pendingNew.a; b = pendingNew.b;
            elseif ~isempty(k) && k <= numel(events)
                a = events(k).start; b = events(k).end;
            else
                set(box,'Visible','off'); continue;
            end
            if isempty(vid), set(box,'Visible','off'); continue; end
            x0 = min(max(a,viewStart),viewEnd); x1 = max(min(b,viewEnd),viewStart);
            set(box,'Position',dataRangeToPix(boxes{q,2},x0,x1),'Visible','on');
        end
        updateSelBar();
    end

    function undo()
        if isempty(undoStack), setStatus('Nothing to undo.','warn'); return; end
        commitEditQuiet();
        events = undoStack{end}; undoStack(end) = [];
        if groupU > 0 && isempty(findUid(groupU)), groupU = 0; end
        currentEvent = []; selStart = NaN; selEnd = NaN;
        eventsChanged(); updateSelectionGraphics();
        setStatus('Undid the last change.','ok');
    end

    % ---------------- playback ----------------------------------------------
    function togglePlay()
        if isPlaying, pausePlay(); else, startPlay(); end
    end

    function setLoop(on)
        loopOn = logical(on);
        for hc = [chkLoop chkBarLoop]
            if ~isempty(hc) && ishghandle(hc), set(hc,'Value',double(loopOn)); end
        end
        if loopOn, setStatus('Loop on: the selection replays until you pause.','info');
        else,      setStatus('Loop off: playback stops at the end of the selection.','info'); end
    end

    function setSpeedIdx(i)
        if isnan(i) || i < 1 || i > numel(SPEEDS), return; end
        playSpeed = SPEEDS(i);
        for hp = [popSpeed popBarSpeed]
            if ~isempty(hp) && ishghandle(hp), set(hp,'Value',i); end
        end
        if isPlaying
            t = NaN; try, t = currentPlayTime(mainPlayer); catch, end
            stopPlay();
            if isfinite(t), cursorTime = clampT(t); end
            startPlay();
        end
        msg = sprintf('Speed %s.',SPEED_NAMES{i});
        if playSpeed < 1, msg = [msg ' Slower speeds also lower the pitch.']; end
        setStatus(msg,'info');
    end

    function onSpeed(src)
        setSpeedIdx(get(src,'Value'));
    end

    function startPlay()
        if isempty(vid) || isPlaying, return; end
        if isempty(mainPlayer), buildPlayer(); end
        if isempty(mainPlayer), return; end
        if hasSel()
            a = selStart; b = selEnd;
            % resume from the play point if it is inside the selection; the
            % selection and the view are never moved
            if ~isnan(cursorTime) && cursorTime > a && cursorTime < b - 0.01
                t0 = cursorTime; whatTxt = 'rest of selection';
            else
                t0 = a; whatTxt = 'selection';
            end
            t1 = b; playLoopA = a;
            buildFrameCache(a,b);        % decode once -> smooth playback and loops
            if a < viewStart || b > viewEnd
                w = max(viewEnd-viewStart, 1.2*(b-a)); c = (a+b)/2; setView(c-w/2,c+w/2);
            end
        elseif ~isnan(cursorTime)
            t0 = cursorTime; t1 = dur; whatTxt = 'from cursor'; playLoopA = NaN;
        else
            t0 = viewStart; t1 = dur; whatTxt = 'from view start'; playLoopA = NaN;
        end
        loopPass = 1;
        s0 = max(1,round(t0*fs)+1); s1 = min(numel(audio),round(t1*fs));
        if s1 <= s0, return; end
        try
            set(mainPlayer,'SampleRate',round(fs*playSpeed));
        catch
            playSpeed = 1; set(popSpeed,'Value',1); set(popBarSpeed,'Value',1);
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
        syncPlayBtn();
        msg = sprintf('Playing %s at %.2fx...   Space = pause.',whatTxt,playSpeed);
        if loopOn && ~isnan(playLoopA), msg = [msg '  Looping (L turns it off).']; end
        if max(abs(audio(s0:s1))) < 1e-4, msg = [msg '  (This part of the recording is silent.)']; end
        setStatus(msg,'info');
    end

    function playTick()
        if ~isPlaying, return; end
        player = mainPlayer;
        if isempty(player), stopPlay(); return; end
        ended = ~isplaying(player);
        if ~ended
            t = currentPlayTime(player);                 % audio = master clock
            ended = t >= playEndTime;
        end
        if ended
            if loopOn && ~isnan(playLoopA), restartLoop(); else, stopPlay(); end
            return;
        end
        if isnan(playLoopA) && (t > viewEnd || t < viewStart)
            % open-ended playback only: page forward a whole window
            w = viewEnd-viewStart; viewStart = max(0,t-0.02*w); viewEnd = min(dur,viewStart+w);
            refreshView();
        end
        if ~isempty(playGfx) && all(ishghandle(playGfx))
            set(playGfx(1),'XData',[t t]); set(playGfx(2),'XData',[t t]);
        end
        if ~isempty(ovCur) && all(isgraphics(ovCur)), set(ovCur,'XData',[t t]); end
        advanceFrameTo(t);
        % plain drawnow: 'drawnow limitrate' caps the screen at 20 fps
        drawnow;
    end

    function stopPlay()
        try, stop(playTimer); end %#ok<TRYNC>
        try, if ~isempty(mainPlayer), stop(mainPlayer); end; catch, end
        isPlaying = false; deleteValid(playGfx); playGfx = gobjects(0);
        if ~isempty(ovCur) && all(isgraphics(ovCur)), set(ovCur,'XData',[cursorTime cursorTime]); end
        syncPlayBtn();
    end

    % ---------------- save / load -------------------------------------------
    function C = buildSaveCells()
        % One row per item, sorted by start-audiofile-relative. Transcript rows
        % are the pieces; fluent words are numbered 1..N (fluent-event-id) and
        % their disfluencies / errors carry the same number.
        hdr = {'start','stop','start-audiofile-relative','stop-audiofile-relative', ...
            'event-type','transcript','fluent-event-id','disfluency','error','notes'};
        if isempty(events), C = hdr; return; end
        [~,ord] = sort([events.start]); ev = events(ord);
        % a transcript piece whose fluent event covers it exactly is saved once,
        % as the fluent row
        isF = strcmp({ev.kind},'fluent'); dup = false(1,numel(ev));
        for r = find(strcmp({ev.kind},'trans'))
            dup(r) = any(isF & [ev.tok] == ev(r).uid & abs([ev.start]-ev(r).start) < 1e-3 & ...
                abs([ev.end]-ev(r).end) < 1e-3);
        end
        ev = ev(~dup);
        fUids = [ev(strcmp({ev.kind},'fluent')).uid];
        C = cell(numel(ev)+1,numel(hdr)); C(1,:) = hdr;
        for r = 1:numel(ev)
            e = ev(r); tr = ''; fid = ''; dc = ''; ec = '';
            switch e.kind
                case 'trans',  ty = 'transcript'; tr = e.transcript;
                case 'fluent', ty = 'fluent'; tr = e.transcript; fid = find(fUids == e.uid,1);
                case 'disfl'
                    ty = 'disfluency'; dc = e.cat;
                    if DISFL_HAS_TRANSCRIPT, tr = e.transcript; end
                otherwise,     ty = 'error'; ec = e.cat;
            end
            if any(strcmp(e.kind,{'disfl','error'})) && e.link > 0
                m = find(fUids == e.link,1); if ~isempty(m), fid = m; end
            end
            C(r+1,:) = {numOrEmpty(e.gStart),numOrEmpty(e.gStop),e.start,e.end, ...
                ty,tr,fid,dc,ec,e.notes};
        end
    end

    function ok = saveAnnotations()
        % .xlsx with three sheets:
        %   events          - what you asked for: start, end, duration, fluent
        %                     event, its disfluency / error codes, notes
        %   codes           - what each code means (R = rep, ...)
        %   annotator_data  - everything the annotator needs to reopen the file
        %                     (transcript pieces included); you can ignore it
        ok = false;
        if isempty(vid), errordlg('Load a video first.','Save'); return; end
        commitEdit();
        if isempty(events)
            c = questdlg('There is nothing to save yet. Save an empty annotation file anyway?', ...
                'Nothing to save','Save anyway','Cancel','Cancel');
            if ~strcmp(c,'Save anyway'), return; end
        end
        [vdir,base] = fileparts(videoPath);
        defName = fullfile(vdir,[base '_annot-disfluencies.xlsx']);
        [fn,fp] = uiputfile({'*.xlsx','Excel workbook (*.xlsx)';'*.csv','CSV files (*.csv)'}, ...
            'Save annotations',defName);
        if isequal(fn,0), return; end
        out = fullfile(fp,fn);
        E = exportCells(); K = codeTable(); A = buildSaveCells();
        [od,ob,oext] = fileparts(out);
        csvMode = strcmpi(oext,'.csv');
        if ~csvMode
            try
                if exist(out,'file'), delete(out); end
                writecell(E,out,'Sheet','events');
                writecell(K,out,'Sheet','codes');
                writecell(A,out,'Sheet','annotator_data');
            catch err
                warndlg(sprintf(['Excel saving failed on this MATLAB (%s), so CSV files are saved ' ...
                    'instead, next to each other in the same folder.'],err.message),'Saved as CSV');
                csvMode = true; out = fullfile(od,[ob '.csv']);
            end
        end
        if csvMode
            [od,ob] = fileparts(out);
            try
                writeCsvManual(out,E(1,:),E(2:end,:));
                writeCsvManual(fullfile(od,[ob '_codes.csv']),K(1,:),K(2:end,:));
                writeCsvManual(fullfile(od,[ob '_annotator.csv']),A(1,:),A(2:end,:));
            catch err
                errordlg(sprintf(['Could not save the annotations:\n\n%s\n\nIf the file is open in ' ...
                    'Excel, close it and try again.'],err.message),'Save error');
                return;
            end
        end
        nF = sum(strcmp({events.kind},'fluent'));
        nD = sum(ismember({events.kind},{'disfl','error'}));
        dirty = false; updateTitle(); ok = true;
        [~,fn1,fe1] = fileparts(out);
        setStatus(sprintf('Saved %d fluent event(s) and %d disfluency / error event(s) to %s%s', ...
            nF,nD,fn1,fe1),'ok');
        if csvMode
            where = sprintf('%s\n(+ %s_codes.csv and %s_annotator.csv)',out,ob,ob);
        else
            where = sprintf('%s\n\nSheets: events, codes (what R / P / ... mean), annotator_data (to reopen it here)',out);
        end
        msgbox(sprintf('Saved %d fluent event(s) and %d disfluency / error event(s) to:\n\n%s', ...
            nF,nD,where),'Annotations saved','help','replace');
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
            C = readSavedTable(f);
        catch err
            errordlg(sprintf(['Could not read "%s".\n\n%s\n\nTry re-saving it from Excel ' ...
                'as .xlsx or .csv and load it again.'],fn,err.message),'Load error');
            return;
        end
        if size(C,1) < 2
            errordlg(sprintf('"%s" has no annotation rows.',fn),'Nothing to load'); return;
        end

        % columns by header name (punctuation / case ignored)
        hdrN = regexprep(lower(cellfun(@toStr,C(1,:),'UniformOutput',false)),'[^a-z0-9]','');
        colOf = @(nm) firstOr0(find(strcmp(hdrN,nm),1));
        cS = colOf('startaudiofilerelative'); cE = colOf('stopaudiofilerelative');
        cTyp = colOf('eventtype'); cTr = colOf('transcript'); cWi = colOf('wordindex');
        cId = colOf('fluenteventid');
        cD = colOf('disfluency'); cEr = colOf('error'); cN = colOf('notes');
        cGS = colOf('start'); cGE = colOf('stop'); cOld = 0;
        if cS > 0 && cE > 0
            if cTyp > 0, fmt = 'typed'; else, fmt = 'v3'; end
        else                                        % v2: starts / ends / event_type
            cS = colOf('starts'); cE = colOf('ends'); cOld = colOf('eventtype');
            cTyp = 0; cGS = 0; cGE = 0; fmt = 'v2';
            if cS == 0 || cE == 0
                errordlg(sprintf(['This file does not look like an annotation file.\n\n' ...
                    'Expected columns: start-audiofile-relative, stop-audiofile-relative, ' ...
                    'event-type, transcript, fluent-event-id, disfluency, error, notes.']),'Invalid file');
                return;
            end
        end
        cellAt = @(i,c) cellOr(C,i,c);

        R = struct('row',{},'st',{},'en',{},'typ',{},'tr',{},'fid',{},'wi',{},'dc',{},'ec',{}, ...
            'nt',{},'gs',{},'ge',{});
        nBad = 0; nOut = 0;
        for i = 2:size(C,1)
            if all(cellfun(@(v)isempty(toStr(v)),C(i,:))), continue; end
            st = toNum(cellAt(i,cS)); en = toNum(cellAt(i,cE));
            if isnan(st) || isnan(en) || en <= st, nBad = nBad + 1; continue; end
            if st >= dur, nOut = nOut + 1; continue; end
            dc = lower(toStr(cellAt(i,cD))); ec = lower(toStr(cellAt(i,cEr)));
            switch fmt
                case 'typed'
                    typ = lower(toStr(cellAt(i,cTyp)));
                case 'v3'
                    if ~isempty(dc) && ~isempty(ec), typ = 'both';
                    elseif ~isempty(dc),             typ = 'disfluency';
                    elseif ~isempty(ec),             typ = 'error';
                    else,                            typ = 'fluent';
                    end
                otherwise
                    typ = 'disfluency';
                    switch lower(toStr(cellAt(i,cOld)))
                        case 'repetition',     dc = 'rep';
                        case 'block',          dc = 'bl';
                        case 'prolongation',   dc = 'prol';
                        case {'text-only',''}, typ = 'transcript';
                        otherwise,             dc = lower(toStr(cellAt(i,cOld)));
                    end
            end
            R(end+1) = struct('row',i,'st',st,'en',min(en,dur),'typ',typ,'tr',toStr(cellAt(i,cTr)), ...
                'fid',toNum(cellAt(i,cId)),'wi',toNum(cellAt(i,cWi)),'dc',dc,'ec',ec, ...
                'nt',toStr(cellAt(i,cN)),'gs',toNum(cellAt(i,cGS)),'ge',toNum(cellAt(i,cGE))); %#ok<AGROW>
        end
        if isempty(R)
            errordlg(sprintf('No usable rows were found in "%s".',fn),'Nothing to load'); return;
        end

        % ---- checks ----
        typs = {R.typ}; rowsN = [R.row]; fidAll = [R.fid];
        hasTr = ~cellfun(@isempty,{R.tr});
        isDi = ismember(typs,{'disfluency','both'}); isEr = ismember(typs,{'error','both'});
        isFl = strcmp(typs,'fluent'); isTx = strcmp(typs,'transcript');
        isTg = strcmp(typs,'target');                % from the short-lived v4 layout
        issues = {};
        if nBad > 0, issues{end+1} = sprintf('- %d row(s) with missing or invalid times were skipped.',nBad); end
        if nOut > 0, issues{end+1} = sprintf('- %d row(s) start after the end of this video and were skipped.',nOut); end
        if strcmp(fmt,'v2')
            issues{end+1} = ['- Old layout (starts / ends / event_type): text-only rows load as ' ...
                'transcript, the rest as unresolved disfluencies.'];
        end
        bad = find(diff([R.st]) < 0);
        if ~isempty(bad)
            issues{end+1} = sprintf('- Rows are not sorted by start-audiofile-relative (rows %s).',rowList(rowsN(bad+1)));
        end
        if strcmp(fmt,'typed')
            b = ~ismember(typs,{'transcript','fluent','disfluency','error','target'});
            if any(b), issues{end+1} = sprintf(['- Unknown event-type (rows %s); expected transcript, ' ...
                    'fluent, disfluency or error. Skipped.'],rowList(rowsN(b))); end
            b = (isTx | isFl) & (~cellfun(@isempty,{R.dc}) | ~cellfun(@isempty,{R.ec}));
            if any(b), issues{end+1} = sprintf(['- Transcript / fluent rows have a disfluency / error ' ...
                    'label (rows %s) - ignored.'],rowList(rowsN(b))); end
        end
        b = isEr & ~isDi & hasTr;
        if any(b), issues{end+1} = sprintf('- Error rows have a transcript (rows %s) - it will be dropped.',rowList(rowsN(b))); end
        if ~DISFL_HAS_TRANSCRIPT
            b = isDi & hasTr;
            if any(b), issues{end+1} = sprintf('- Disfluency rows have a transcript (rows %s).',rowList(rowsN(b))); end
        end
        b = isDi & ~ismember({R.dc},{DISFL.abbr});
        if any(b), issues{end+1} = sprintf('- Unknown disfluency labels (rows %s); expected %s.', ...
                rowList(rowsN(b)),strjoin({DISFL.abbr},', ')); end
        b = isEr & ~ismember({R.ec},{ERRS.abbr});
        if any(b), issues{end+1} = sprintf('- Unknown error labels (rows %s); expected %s.', ...
                rowList(rowsN(b)),strjoin({ERRS.abbr},', ')); end
        b = isFl & ~hasTr;
        if any(b), issues{end+1} = sprintf('- Fluent rows with no word (rows %s).',rowList(rowsN(b))); end
        flIds = fidAll(isFl);
        if any(strcmp(fmt,{'typed','v3'}))
            if any(isnan(flIds))
                issues{end+1} = sprintf('- Fluent rows without a fluent-event-id (rows %s).', ...
                    rowList(rowsN(isFl & isnan(fidAll))));
            end
            v = flIds(~isnan(flIds));
            if numel(unique(v)) < numel(v)
                issues{end+1} = '- Fluent-event-ids are repeated across fluent rows.';
            elseif ~isequal(flIds,1:numel(flIds))
                issues{end+1} = '- Fluent-event-ids are not sequential integers from 1 (renumbered on save).';
            end
            b = (isDi | isEr) & ~isnan(fidAll) & ~ismember(fidAll,flIds);
            if any(b), issues{end+1} = sprintf(['- Disfluency / error rows point to a fluent-event-id ' ...
                    'that no fluent row has (rows %s) - loaded unresolved.'],rowList(rowsN(b))); end
        end
        txSpans = [[R(isTx).st]; [R(isTx).en]];
        inTx = false(1,numel(R));
        for r = find(isFl)
            inTx(r) = any(txSpans(1,:) <= R(r).st + 1e-6 & txSpans(2,:) >= R(r).en - 1e-6);
        end
        if ~isempty(issues)
            c = questdlg(sprintf('Checks found problems in "%s":\n\n%s\n\nLoad it anyway?', ...
                fn,strjoin(issues,sprintf('\n'))),'Annotation file check','Load anyway','Cancel','Cancel');
            if ~strcmp(c,'Load anyway'), setStatus('Load cancelled.','info'); return; end
        end

        % ---- build ----
        newEv = EV0;
        for r = find(isTx)
            newEv(end+1) = mkEv(R(r),'trans','',R(r).tr,0); %#ok<AGROW>
        end
        for r = find(isTg)                           % v4 target words: keep their text as transcript
            if ~any(txSpans(1,:) < R(r).en & txSpans(2,:) > R(r).st)
                newEv(end+1) = mkEv(R(r),'trans','',R(r).tr,0); %#ok<AGROW>
            end
        end
        idKeys = []; idUids = [];
        for r = find(isFl)
            if ~inTx(r)                              % fluent word needs a transcript around it
                newEv(end+1) = mkEv(R(r),'trans','',R(r).tr,0); %#ok<AGROW>
            end
            newEv(end+1) = mkEv(R(r),'fluent','',normWordKeep(R(r).tr),0); %#ok<AGROW>
            newEv(end).wi = pickWordIdx(numel(newEv),R(r).wi);
            if ~isnan(R(r).fid) && ~ismember(R(r).fid,idKeys)
                idKeys(end+1) = R(r).fid; idUids(end+1) = newEv(end).uid; %#ok<AGROW>
            end
        end
        for r = find(isDi | isEr)
            lk = 0;
            if ~strcmp(fmt,'v2') && ~isnan(R(r).fid)
                m = find(idKeys == R(r).fid,1); if ~isempty(m), lk = idUids(m); end
            end
            if isDi(r)
                tr = ''; if DISFL_HAS_TRANSCRIPT, tr = R(r).tr; end
                newEv(end+1) = mkEv(R(r),'disfl',R(r).dc,tr,lk); %#ok<AGROW>
            end
            if isEr(r)
                newEv(end+1) = mkEv(R(r),'error',R(r).ec,'',lk); %#ok<AGROW>
            end
        end

        % events that sit exactly on a transcript piece belong to it
        tr = find(strcmp({newEv.kind},'trans'));
        for j = find(~strcmp({newEv.kind},'trans'))
            m = tr(abs([newEv(tr).start]-newEv(j).start) < 1e-3 & abs([newEv(tr).end]-newEv(j).end) < 1e-3);
            if ~isempty(m), newEv(j).tok = newEv(m(1)).uid; end
        end

        % ---- merge or replace ----
        mergeMode = false;
        if ~isempty(events)
            c = questdlg(sprintf(['"%s" has %d item(s). This video already has %d.\n\n' ...
                'Merge them into your current annotations, or replace the current ones? ' ...
                'Ctrl+Z undoes either.'],fn,numel(newEv),numel(events)),'Load annotations', ...
                'Merge into current','Replace current','Cancel','Merge into current');
            switch c
                case 'Merge into current'
                    mergeMode = true;
                    [isDup,dupUid] = compareToCurrent(newEv);
                    for i = find(isDup)            % re-point links to the existing copy
                        for j = 1:numel(newEv)
                            if newEv(j).link == newEv(i).uid, newEv(j).link = dupUid(i); end
                        end
                    end
                    newEv = newEv(~isDup);
                    if isempty(newEv)
                        msgbox('Everything in that file is already on the timeline - nothing was added.', ...
                            'Nothing new','help','replace');
                        return;
                    end
                case 'Replace current'
                otherwise
                    return;
            end
        end

        commitEditQuiet(); pushUndo(); groupU = 0;
        if mergeMode
            events = [events newEv]; dirty = true; verb = 'Merged';
        else
            events = newEv; dirty = false; verb = 'Loaded';
        end
        currentEvent = []; selStart = NaN; selEnd = NaN;
        updateEventList(); updateTitle(); refreshView();
        setStatus(sprintf('%s %d item(s) from %s.',verb,numel(newEv),fn),'ok');

        function ev = mkEv(rr,kind,cat,tr,link)
            ev = struct('start',rr.st,'end',rr.en,'kind',kind,'cat',cat,'transcript',tr, ...
                'notes',rr.nt,'uid',nextUid,'link',link,'wi',0,'tok',0,'gStart',rr.gs,'gStop',rr.ge);
            nextUid = nextUid + 1;
        end

        function wi = pickWordIdx(fk,given)
            % word-index from the file if it matches, else the first matching
            % word of the surrounding transcript that is not already taken
            wi = 0;
            st = newEv(fk).start; en = newEv(fk).end;
            sk2 = find(strcmp({newEv.kind},'trans') & [newEv.start] <= st + 1e-6 & ...
                [newEv.end] >= en - 1e-6,1);
            if isempty(sk2), return; end
            toks = cellfun(@normWord,regexp(newEv(sk2).transcript,'\S+','match'),'UniformOutput',false);
            w = normWord(newEv(fk).transcript);
            if ~isnan(given) && given >= 1 && given <= numel(toks) && strcmp(toks{given},w)
                wi = given; return;
            end
            taken = [];
            for g = find(strcmp({newEv.kind},'fluent'))
                if g ~= fk && newEv(g).start >= newEv(sk2).start - 1e-6 && ...
                        newEv(g).end <= newEv(sk2).end + 1e-6
                    taken(end+1) = newEv(g).wi; %#ok<AGROW>
                end
            end
            m = setdiff(find(strcmp(toks,w)),taken);
            if ~isempty(m), wi = m(1); end
        end
    end

    % ---------------- help windows ------------------------------------------
    function showHelp()
        msg = sprintf([ ...
            'THE WORKFLOW  (example "the t-t-table is on the ground")\n' ...
            '  1. TRANSCRIPT: highlight the speech on the spectrogram, DOUBLE-CLICK the\n' ...
            '     highlight, type   the t-t-table is on the ground   and press Enter.\n' ...
            '     It is split into pieces on the Transcript strip (drag the edges\n' ...
            '     between pieces to fit). This works any time, even in step 3.\n' ...
            '  2. FLUENT WORD: highlight "table" (or click its piece) and press Enter.\n' ...
            '     The word is marked "fluent event selected".\n' ...
            '  3. DISFLUENT PARTS: highlight "t-t" (or one "t-", or a syllable inside a\n' ...
            '     word) and press Enter. A list opens - tick one or more:\n' ...
            '       R rep   P prol   B bl   T td   /   1-6 the error types\n' ...
            '     then Enter. Repeat for each part. W = the whole episode (from the\n' ...
            '     first disfluent part to the end of the word) - e.g. a block over it.\n' ...
            '  4. Esc = done with this word. Next word: back to step 2.\n' ...
            '  Shortcuts: R / P / B / T / E add a type straight away (no list). Enter on\n' ...
            '  a selected disfluency changes its types. Outside a word, R / P / B / T\n' ...
            '  give an unresolved disfluency.\n\n' ...
            'WHAT YOU SEE\n' ...
            '  The Disfl/Error strip shows only the type (rep, td, er-s ...), at its time.\n' ...
            '  The word piece shows a colour band for every type tied to it.\n' ...
            '  Long text shows its first word + "..."; click the item to read it all in\n' ...
            '  the panel left of the video.\n\n' ...
            'NOTES\n' ...
            '  Double-click a highlight on the Notes strip, or click an item and press N\n' ...
            '  (or double-click it on the Notes strip), then type. The note sits at that\n' ...
            '  item''s time. The typing box is wide so you can see it all.\n\n' ...
            'PLAYBACK\n' ...
            '  Space plays the highlight; pause keeps your place. L loops, 1-4 = speed.\n' ...
            '  The highlighted video is loaded into memory first, so it plays smoothly.\n' ...
            '  Scroll on the video to zoom, drag to move, double-click / Fit to reset.\n' ...
            '  "Pop out video" puts it in its own resizable window.\n\n' ...
            'EDITING\n' ...
            '  Double-click a piece to edit its text. Delete removes the selected item.\n' ...
            '  Right-click (two-finger click) for a menu. Ctrl+Z undoes. IPA buttons\n' ...
            '  insert symbols while typing.\n\n' ...
            'EVENTS TABLE (right)\n' ...
            '  One row per fluent event (start, end, word, notes). Click it: the word is\n' ...
            '  highlighted in the transcript and its disfluency / error events open\n' ...
            '  underneath, colour coded (click again to close). Click one of those to\n' ...
            '  select it. Unresolved events are listed at the bottom.\n\n' ...
            'FILES\n' ...
            '  Ctrl+S saves an Excel file with three sheets:\n' ...
            '   events: start, end, duration, fluent event, event code, notes - each\n' ...
            '           fluent event (code F) followed by its disfluency / error rows\n' ...
            '   codes:  what each code means (R rep, P prol, B bl, T td, E-OS ...)\n' ...
            '   annotator_data: used to reopen the file here (includes the transcript)\n' ...
            '  Loading checks the table first.']);
        msgbox(msg,[APPNAME ' - Help'],'help','replace');
    end

    function showShortcuts()
        old = findall(0,'Tag','DA_Shortcuts');
        if ~isempty(old), figure(old(1)); return; end
        rows = { ...
            'H','ANNOTATING','',[]; ...
            'K','Highlight new speech + double-click','Type the transcript (t-t-table ...)',[]; ...
            'K','Double-click highlight on Notes','Type notes for it',[]; ...
            'K','Highlight a word + Enter','Fluent event selected',[]; ...
            'K','Then highlight a part + Enter','Tick type(s) in the list, Enter',[]; ...
            'K','In the list: R P B T / 1-6','Tick disfluency / error types',[]; ...
            'K','W','Whole episode (first part to end of word)',[]; ...
            'K','Esc','Done with this word',[]; ...
            'K','R / P / B / T / E','Add a type straight away',[]; ...
            'K','N','Notes for the selected item',[]; ...
            'K','Double-click an item','Edit text',[]; ...
            'K','Right-click / 2-finger click','Menu',[]; ...
            'H','COLOURS','',[]; ...
            'C','R','rep - repetition (SLD)',DISFL(1).color; ...
            'C','P','prol - prolongation (SLD)',DISFL(2).color; ...
            'C','B','bl - block (SLD)',DISFL(3).color; ...
            'C','T','td - typical disfluency',DISFL(4).color; ...
            'C','E','error (all types)',ERR_COLOR; ...
            'H','PLAYBACK','',[]; ...
            'K','Space','Play / pause the highlight',[]; ...
            'K','L','Loop on / off',[]; ...
            'K','1 / 2 / 3 / 4','Speed 1x / 0.75x / 0.5x / 0.25x',[]; ...
            'K','Scroll on the video','Zoom (drag = move, double-click = fit)',[]; ...
            'H','EDITING & NAVIGATION','',[]; ...
            'K','Delete','Delete the selected item',[]; ...
            'K','Ctrl/Cmd + Z / S','Undo / save',[]; ...
            'K','Left / Right','Step one video frame (Shift: 1 s)',[]; ...
            'K','Mouse wheel','Scroll (Ctrl/Cmd + wheel zooms)',[]; ...
            'K','Ctrl/Cmd + I / O / N / A','Zoom in / out / to selection / all',[]; ...
            'K','F1','Full help',[]};
        lineH = 21; gapH = 8; sw = 640; colKey = 30; colDesc = 320;
        totalH = 16;
        for q = 1:size(rows,1)
            if strcmp(rows{q,1},'H') && q > 1, totalH = totalH + gapH; end
            totalH = totalH + lineH;
        end
        btnH = 54; sh = totalH + btnH + 10;
        mp = getpixelposition(fig);
        sf = figure('Name','Keyboard & mouse shortcuts','NumberTitle','off', ...
            'MenuBar','none','ToolBar','none','Resize','off','Color',[1 1 1], ...
            'Units','pixels','Position',[mp(1)+(mp(3)-sw)/2 mp(2)+(mp(4)-sh)/2 sw sh], ...
            'Tag','DA_Shortcuts','HandleVisibility','off','KeyPressFcn',@closeOnKey);
        ax = axes('Parent',sf,'Units','pixels','Position',[0 btnH sw totalH+10], ...
            'XLim',[0 sw],'YLim',[0 totalH+10],'YDir','reverse','Visible','off', ...
            'HandleVisibility','off');
        y = 14;
        for q = 1:size(rows,1)
            switch rows{q,1}
                case 'H'
                    if q > 1, y = y + gapH; end
                    text(ax,24,y,rows{q,2},'FontWeight','bold','FontSize',11, ...
                        'Color',[0 0 0],'VerticalAlignment','top');
                    line(ax,[24 sw-24],[y+lineH-3 y+lineH-3],'Color',[0.8 0.8 0.8]);
                case 'K'
                    text(ax,colKey,y,rows{q,2},'FontWeight','bold','FontSize',10, ...
                        'Color',[0 0 0],'VerticalAlignment','top');
                    text(ax,colDesc,y,rows{q,3},'FontSize',10,'Color',[0 0 0], ...
                        'VerticalAlignment','top');
                case 'C'
                    patch(ax,colKey+[0 14 14 0],y+[3 3 15 15],rows{q,4}, ...
                        'EdgeColor',rows{q,4}*0.8,'FaceAlpha',0.7);
                    text(ax,colKey+22,y,rows{q,2},'FontWeight','bold','FontSize',10, ...
                        'Color',[0 0 0],'VerticalAlignment','top');
                    text(ax,colDesc,y,rows{q,3},'FontSize',10,'Color',[0 0 0], ...
                        'VerticalAlignment','top');
            end
            y = y + lineH;
        end
        uicontrol(sf,'Style','pushbutton','String','Close','Units','pixels', ...
            'Position',[(sw-110)/2 14 110 30],'Callback',@(~,~)delete(sf));

        function closeOnKey(src,e)
            if any(strcmp(e.Key,{'escape','return','space'})), delete(src); end
        end
    end


    function s = deLabel(k)
        % The Disfl/Error strip shows only the type, e.g.  rep  /  er-s
        s = events(k).cat;
    end

    function computeLanes()
        % Overlapping items stack in rows (lanes) on the Disfl/Error strip:
        % fluent / disfluent pieces in the upper part, errors in the lower part.
        laneOf = zeros(1,numel(events)); nLaneP = 1; nLaneE = 1;
        if isempty(events), return; end
        P = find(strcmp({events.kind},'disfl'));
        [l,n] = assignLanes(P); laneOf(P) = l; nLaneP = max(1,n);
        E = find(strcmp({events.kind},'error'));
        [l,n] = assignLanes(E); laneOf(E) = l; nLaneE = max(1,n);
    end

    function [lane,n] = assignLanes(idx)
        lane = zeros(1,numel(idx)); n = 0;
        if isempty(idx), return; end
        st = [events(idx).start]; du = [events(idx).end] - st;
        [~,o] = sortrows([st(:) -du(:)]);          % earlier first; longer first on ties
        ends = [];
        for j = reshape(o,1,[])
            L = find(ends <= events(idx(j)).start + TOL,1);
            if isempty(L), ends(end+1) = events(idx(j)).end; L = numel(ends); %#ok<AGROW>
            else, ends(L) = events(idx(j)).end; end
            lane(j) = L;
        end
        n = numel(ends);
    end

    function bnd = laneBand(k)
        if numel(laneOf) ~= numel(events), computeLanes(); end
        if strcmp(events(k).kind,'error'), r = [0.02 0.34]; n = nLaneE;
        else,                              r = [0.36 0.98]; n = nLaneP; end
        L = max(1,laneOf(k)); h = (r(2)-r(1))/n;
        top = r(2) - (L-1)*h;
        bnd = [top-h+0.008, top-0.008];
    end

    function onMultiKey()
        openEventPicker();
    end

    function k = sameWindow(a,b,kind,cat)
        k = [];
        if isempty(events), return; end
        k = find(strcmp({events.kind},kind) & strcmp({events.cat},cat) & ...
            abs([events.start]-a) < 1e-3 & abs([events.end]-b) < 1e-3,1);
    end

    function u = linkFor(k)
        % The fluent word a disfluency / error belongs to: the one it overlaps,
        % else the next fluent word starting within 1.5 s after it.
        u = 0;
        fl = find(strcmp({events.kind},'fluent'));
        if isempty(fl), return; end
        a = events(k).start; b = events(k).end;
        ov = min(b,[events(fl).end]) - max(a,[events(fl).start]);
        if any(ov > TOL)
            [~,m] = max(ov); u = events(fl(m)).uid; return;
        end
        after = fl([events(fl).start] >= b - TOL & [events(fl).start] - b <= 1.5);
        if ~isempty(after), [~,m] = min([events(after).start]); u = events(after(m)).uid; end
    end

    function relinkAll()
        % Keep existing connections; connect the rest to the fluent word they
        % overlap, or the next one after them.
        for j = find(ismember({events.kind},{'disfl','error'}))
            if events(j).link > 0 && ~isempty(findUid(events(j).link)), continue; end
            events(j).link = linkFor(j);
        end
    end


    function pc = pieceColors(k)
        % Colour bands for transcript piece k: disfluency / error marks on it
        % (or over it), plus - for a fluent word - everything that resolves to
        % it. Fluent itself gets no band - it is implied.
        pc = zeros(0,3); keys = {};
        own = attachedTo(k);
        a = events(k).start; b = events(k).end;
        over = find(ismember({events.kind},{'disfl','error'}) & [events.tok] == 0 & ...
            (min(b,[events.end]) - max(a,[events.start])) >= 0.5*(b-a));
        m = [own over];
        for f = own(strcmp({events(own).kind},'fluent'))
            ln = linkedAll(events(f).uid);
            [~,o] = sort([events(ln).start]);
            m = [m ln(o)]; %#ok<AGROW>
        end
        for j = m
            if strcmp(events(j).kind,'fluent'), continue; end
            if strcmp(events(j).kind,'error'), key = 'error'; else, key = [events(j).kind events(j).cat]; end
            if any(strcmp(keys,key)), continue; end
            keys{end+1} = key; pc(end+1,:) = eventColor(j); %#ok<AGROW>
        end
    end

    function idx = attachedTo(k)
        % Events sitting on transcript piece k.
        idx = [];
        if isempty(k) || k > numel(events) || ~strcmp(events(k).kind,'trans'), return; end
        idx = find([events.tok] == events(k).uid & ~strcmp({events.kind},'trans'));
    end

    function syncTok(k)
        % Events on a transcript piece follow it when it is resized / moved.
        for j = attachedTo(k)
            events(j).start = events(k).start; events(j).end = events(k).end; clearGlobal(j);
        end
    end

    function rightClickPick(t,panel,y)
        % Right-click target: the highlight if you clicked inside it, else the
        % item under the pointer (a transcript piece on the Transcript strip).
        if hasSel() && isempty(currentEvent) && t >= selStart && t <= selEnd
            return;                                   % right-click inside a highlight
        end
        idx = eventAtTime(t,panel,y);
        if ~isempty(idx) && ~isequal(idx,currentEvent)
            currentEvent = idx; selStart = events(idx).start; selEnd = events(idx).end;
            cursorTime = selStart; syncListSelection();
            redrawEvents(); updateSelectionGraphics(); updateCursorGraphics();
        end
    end

    function pickFromPointer()
        % Same as above, worked out from where the pointer is (in case the menu
        % opens before the click handler has run).
        axs = [axTrans axDE axNotes axSpec axWave];
        pnl = {'transcript','disferr','notes','',''};
        for q = 1:numel(axs)
            cp = get(axs(q),'CurrentPoint'); xl = get(axs(q),'XLim'); yl = get(axs(q),'YLim');
            if cp(1,1) >= xl(1) && cp(1,1) <= xl(2) && cp(1,2) >= yl(1) && cp(1,2) <= yl(2)
                rightClickPick(clampT(cp(1,1)),pnl{q},cp(1,2)); return;
            end
        end
    end

    function j = touchingPiece(k,side)
        % The transcript piece that shares the edge being dragged.
        j = [];
        tr = find(strcmp({events.kind},'trans')); tr = tr(tr ~= k);
        if isempty(tr), return; end
        if strcmp(side,'left'), j = tr(find(abs([events(tr).end] - events(k).start) < 1e-3,1));
        else,                   j = tr(find(abs([events(tr).start] - events(k).end) < 1e-3,1)); end
    end

    function [a,b,tk] = currentTarget()
        % What the menu / keys act on: the selected transcript piece, the piece
        % a selected event sits on, the selected event's window, or the highlight.
        a = NaN; b = NaN; tk = [];
        if ~isempty(currentEvent) && currentEvent <= numel(events)
            ev = events(currentEvent);
            if strcmp(ev.kind,'trans'), a = ev.start; b = ev.end; tk = currentEvent; return; end
            if ev.tok > 0
                t2 = findUid(ev.tok);
                if ~isempty(t2), a = events(t2).start; b = events(t2).end; tk = t2; return; end
            end
            a = ev.start; b = ev.end; return;
        end
        if hasSel()
            a = selStart; b = selEnd;
            tk = find(strcmp({events.kind},'trans') & abs([events.start]-a) < 1e-3 & ...
                abs([events.end]-b) < 1e-3,1);
        end
    end

    function j = hasAnnot(a,b,kind,cat)
        j = [];
        if isempty(events), return; end
        j = find(strcmp({events.kind},kind) & strcmp({events.cat},cat) & ...
            abs([events.start]-a) < 1e-3 & abs([events.end]-b) < 1e-3,1);
    end

    function onMenuOpen(menu) %#ok<INUSD>
        if isEditing, commitEdit(); end
        pickFromPointer();
        [a,b,tk] = currentTarget();
        kids = [miFlu miDis miErrRoot miSev miSplit miEdit miPlay miDel];
        if isnan(a)
            set(kids,'Enable','off');
            setStatus('Highlight some speech or click a transcript piece first, then right-click.','warn');
            return;
        end
        set(kids,'Enable','on');
        onoff = {'off','on'};
        set(miFlu,'Checked',onoff{1+~isempty(hasAnnot(a,b,'fluent',''))});
        for q = 1:numel(DISFL)
            set(miDis(q),'Checked',onoff{1+~isempty(hasAnnot(a,b,'disfl',DISFL(q).abbr))});
        end
        anyErr = false;
        for q = 1:numel(ERRS)
            on = ~isempty(hasAnnot(a,b,'error',ERRS(q).abbr)); anyErr = anyErr || on;
            set(miErr(q),'Checked',onoff{1+on});
        end
        try, set(miErrRoot,'Checked',onoff{1+anyErr}); catch, end
        canSplit = ~isempty(tk) && numel(splitTokens(events(tk).transcript)) > 1;
        set(miSplit,'Enable',onoff{1+canSplit});
        set([miEdit miDel],'Enable',onoff{1+~isempty(currentEvent)});
        if ~isempty(tk)
            set(miEdit,'Label',sprintf('Edit text  ("%s")',events(tk).transcript));
        else
            set(miEdit,'Label','Edit text');
        end
    end

    function menuToggle(kind,cat)
        % Tick = add that annotation to the piece; untick = remove it.
        [a,b,tk] = currentTarget();
        if isnan(a), return; end
        j = hasAnnot(a,b,kind,cat);
        if isempty(j)
            addAnnotation(kind,cat);
        else
            keepU = 0; if ~isempty(tk), keepU = events(tk).uid; end
            pushUndo(); events(j) = []; relinkAll();
            if keepU > 0, currentEvent = findUid(keepU); else, currentEvent = []; end
            eventsChanged(); updateSelectionGraphics();
            setStatus(sprintf('Removed it. This piece now has: %s.',pieceTypes(a,b)),'ok');
        end
    end

    function k = addAnnotation(kind,cat)
        % Add one annotation on the current target (piece / event / highlight).
        k = [];
        [a,b,tk] = currentTarget();
        if isnan(a)
            setStatus('Highlight some speech or click a transcript piece first.','warn'); return;
        end
        if ~isempty(hasAnnot(a,b,kind,cat))
            setStatus(sprintf('This piece already has that. It has: %s.',pieceTypes(a,b)),'info'); return;
        end
        if strcmp(kind,'fluent')
            ov = overlapIdx(a,b,{'fluent'},[]);
            if ~isempty(ov)
                setStatus(sprintf('That overlaps the fluent word "%s".',events(ov(1)).transcript),'warn');
                beep; return;
            end
        end
        pushUndo();
        k = addRaw(a,b,tk,kind,cat);
        relinkAll();
        if ~isempty(tk), currentEvent = tk; else, currentEvent = k; end
        selStart = a; selEnd = b;
        eventsChanged(); updateSelectionGraphics();
        setStatus(sprintf('This piece has: %s. Right-click again to add more (Ctrl+Z to undo).', ...
            pieceTypes(a,b)),'ok');
        if isempty(tk) && DISFL_HAS_TRANSCRIPT && ~strcmp(kind,'error') && ~strcmp(cat,'bl')
            % a bare highlight has no text yet: type it
            startEditing(k,'text');
            if strcmp(kind,'fluent'), what = 'the word (e.g. table)'; else, what = 'what was said (e.g. t)'; end
            setStatus(sprintf('Added %s - type %s, Enter saves, Esc skips.',pieceTypes(a,b),what),'info');
        end
    end

    function k = addRaw(a,b,tk,kind,cat)
        txt = '';
        if ~isempty(tk) && ~strcmp(kind,'error'), txt = tokenWord(events(tk).transcript); end
        lk = 0;
        if groupU > 0 && any(strcmp(kind,{'disfl','error'})), lk = groupU; end
        k = newEvent(a,b,kind,cat,txt,lk);
        if ~isempty(tk), events(k).tok = events(tk).uid; end
    end

    function s = pieceTypes(a,b)
        j = find(~strcmp({events.kind},'trans') & abs([events.start]-a) < 1e-3 & abs([events.end]-b) < 1e-3);
        if isempty(j), s = 'nothing'; return; end
        labs = cell(1,numel(j));
        for q = 1:numel(j)
            if strcmp(events(j(q)).kind,'fluent'), labs{q} = 'fluent'; else, labs{q} = events(j(q)).cat; end
        end
        s = strjoin(labs,' + ');
    end

    function w = tokenWord(s)
        % "t-" -> "t", "table," -> "table"
        w = regexprep(strtrim(char(s)),'^[-\.,;:!?"''()]+|[-\.,;:!?"''()]+$','');
        if isempty(w), w = strtrim(char(s)); end
    end

    function toks = splitTokens(str)
        % "t-t-table is" -> {'t-','t-','table','is'}
        toks = regexp(char(str),'[^\s-]+-*','match');
        if isempty(toks), toks = {strtrim(char(str))}; end
    end

    function k1 = addPieces(a,b,str)
        % One transcript piece per word / part-word, spread over [a b] in
        % proportion to its length (drag the edges to fit the audio).
        toks = splitTokens(str);
        w = max(1,cellfun(@numel,toks));
        edges = a + (b-a)*[0 cumsum(w)]/sum(w);
        k1 = [];
        for i = 1:numel(toks)
            k = newEvent(edges(i),edges(i+1),'trans','',toks{i},0);
            if isempty(k1), k1 = k; end
        end
    end

    function menuSplit()
        [~,~,tk] = currentTarget();
        if isempty(tk), return; end
        ev = events(tk);
        if numel(splitTokens(ev.transcript)) < 2, return; end
        pushUndo();
        for j = attachedTo(tk), events(j).tok = 0; end     % they keep their window
        events(tk) = [];
        k1 = addPieces(ev.start,ev.end,ev.transcript);
        currentEvent = k1; selStart = events(k1).start; selEnd = events(k1).end;
        eventsChanged(); updateSelectionGraphics();
        setStatus('Split into pieces - drag the edges between them to fit the audio.','ok');
    end

    function startPlayFresh()
        stopPlay();
        [a,b] = currentTarget();
        if ~isnan(a), selStart = a; selEnd = b; end
        cursorTime = selStart; updateSelectionGraphics(); startPlay();
    end

    function resetVideoZoom()
        vzoom = []; applyVideoZoom();
    end

    % ---------------- video display, zoom and smooth playback ----------------
    function r = videoRect()
        % Part of the frame on show, in frame pixels [x0 x1 y0 y1].
        if isempty(vzoom), r = [1 vid.Width 1 vid.Height]; else, r = vzoom; end
    end

    function updateFrameStride()
        % Show a decimated frame when the shown region is much larger than the
        % panel: far less pixel data to push to the screen on every tick.
        frameStride = 1;
        if isempty(vid), return; end
        p = getpixelposition(axVideo,true); r = videoRect();
        frameStride = max(1, floor(min((r(4)-r(3)+1)/max(p(4),1), (r(2)-r(1)+1)/max(p(3),1))));
    end

    function fr = dispFrame(fr)
        r = videoRect();
        r(2) = min(r(2),size(fr,2)); r(4) = min(r(4),size(fr,1));
        fr = fr(r(3):frameStride:r(4), r(1):frameStride:r(2), :);
    end

    function setFrame(fullFr)
        lastFull = fullFr;
        setImgData(dispFrame(fullFr));
    end

    function setImgData(d)
        if isempty(imgVideo) || ~ishghandle(imgVideo), return; end
        r = videoRect();
        set(imgVideo,'CData',d,'XData',[r(1) r(2)],'YData',[r(3) r(4)]);
    end

    function advanceFrameTo(t)
        % Playback frame update: from the in-memory cache when the played
        % stretch is cached (smooth), else decode forward from the file.
        if isempty(vid) || isempty(imgVideo) || ~ishghandle(imgVideo), return; end
        if cacheCovers(t)
            i = find(fcTimes <= t + 1e-6,1,'last'); if isempty(i), i = 1; end
            if i ~= fcLast, setImgData(fcFrames(:,:,:,i)); fcLast = i; end
            return;
        end
        fp = 1/max(vid.FrameRate,1);
        try
            if t < vid.CurrentTime - fp || t > vid.CurrentTime + 0.5
                vid.CurrentTime = max(0, min(t, vid.Duration - fp));
            end
            fr = [];
            while hasFrame(vid) && vid.CurrentTime <= t
                fr = readFrame(vid);
            end
            if ~isempty(fr), setFrame(fr); end
        catch
        end
    end

    function showFrameAt(t)
        % Random-access frame display (clicks, drags, pauses).
        if isempty(vid) || isempty(imgVideo) || ~ishghandle(imgVideo) || isnan(t), return; end
        if cacheCovers(t)
            i = find(fcTimes <= t + 1e-6,1,'last'); if isempty(i), i = 1; end
            setImgData(fcFrames(:,:,:,i)); fcLast = i; return;
        end
        tt = max(0,min(t, vid.Duration - 1/max(vid.FrameRate,1)));
        try
            vid.CurrentTime = tt; setFrame(readFrame(vid));
        catch
        end
    end

    function tf = cacheCovers(t)
        tf = ~isempty(fcTimes) && isequal(fcZoom,vzoom) && fcStride == frameStride && ...
            t >= fcTimes(1) - 1e-6 && t <= fcTimes(end) + 1/max(fcFps,1);
    end

    function clearFrameCache()
        fcFrames = []; fcTimes = []; fcRange = []; fcLast = 0;
    end

    function ok = buildFrameCache(a,b)
        % Decode the frames for [a b] once, into memory, so playing / looping
        % that stretch never has to decode on the fly.
        ok = false;
        if isempty(vid), return; end
        fps = max(vid.FrameRate,1);
        if ~isempty(fcTimes) && fcRange(1) <= a + 1e-6 && fcRange(2) >= b - 1e-6 && ...
                isequal(fcZoom,vzoom) && fcStride == frameStride
            ok = true; return;
        end
        if b - a > CACHE_MAX_S, return; end
        r = videoRect();
        h = numel(r(3):frameStride:r(4)); w = numel(r(1):frameStride:r(2));
        n = ceil((b-a)*fps) + 3;
        if n*h*w*3/1e6 > CACHE_MAX_MB, return; end
        clearFrameCache();
        setStatus(sprintf('Preparing %.1f s of video for smooth playback...',b-a),'info'); drawnow;
        frames = []; times = zeros(1,n); c = 0;
        try
            vid.CurrentTime = max(0, a - 1/fps);
            while hasFrame(vid) && c < n
                tt = vid.CurrentTime; fr = readFrame(vid);
                if tt + 1/fps < a, continue; end
                d = dispFrame(fr);
                if isempty(frames), frames = zeros([size(d) n],'like',d); end
                c = c + 1; frames(:,:,:,c) = d; times(c) = tt;
                lastFull = fr;
                if tt > b, break; end
            end
        catch
        end
        if c == 0, return; end
        fcFrames = frames(:,:,:,1:c); fcTimes = times(1:c); fcRange = [a b];
        fcZoom = vzoom; fcStride = frameStride; fcFps = fps; fcLast = 0;
        ok = true;
    end

    function tf = overVideo()
        % pointer over the docked video (the popped-out window has its own scroll)
        tf = false;
        if isempty(axVideo) || ~ishghandle(axVideo) || ~isequal(ancestor(axVideo,'figure'),fig), return; end
        p = getpixelposition(axVideo,true); fp = get(fig,'CurrentPoint');
        tf = fp(1) >= p(1) && fp(1) <= p(1)+p(3) && fp(2) >= p(2) && fp(2) <= p(2)+p(4);
    end

    function videoZoomBy(f,cx,cy)
        % Zoom the video picture by factor f around frame point (cx, cy).
        if isempty(vid) || isempty(imgVideo), return; end
        W = vid.Width; H = vid.Height; r = videoRect();
        if nargin < 3, cx = (r(1)+r(2))/2; cy = (r(3)+r(4))/2; end
        w = (r(2)-r(1)+1)/f;
        if w >= W - 1
            vzoom = [];
        else
            w = max(round(w),40); h = max(round(w*H/W),20);
            x0 = round(min(max(cx - w/2,1), W - w + 1));
            y0 = round(min(max(cy - h/2,1), H - h + 1));
            vzoom = [x0 x0+w-1 y0 y0+h-1];
        end
        applyVideoZoom();
    end

    function applyVideoZoom()
        clearFrameCache(); updateFrameStride();
        r = videoRect();
        set(axVideo,'XLim',[r(1)-0.5 r(2)+0.5],'YLim',[r(3)-0.5 r(4)+0.5]);
        if ~isempty(lastFull), setImgData(dispFrame(lastFull)); end
        if isempty(vzoom)
            setStatus('Video zoom reset (whole picture).','info');
        else
            setStatus(sprintf(['Video zoom %.1fx - scroll on the video to zoom, drag it to move ' ...
                'around, double-click or Fit to reset.'],vid.Width/(r(2)-r(1)+1)),'info');
        end
    end

    function onVideoDown()
        if isempty(vid) || isempty(imgVideo), return; end
        vf = ancestor(axVideo,'figure');
        if strcmp(get(vf,'SelectionType'),'open'), vzoom = []; applyVideoZoom(); return; end
        if isempty(vzoom)
            setStatus('Scroll on the video (or click +) to zoom in; then drag to move around.','info');
            return;
        end
        clearFrameCache();
        vpanStart = get(vf,'CurrentPoint'); vpanZoom = vzoom;
        set(vf,'WindowButtonMotionFcn',@(~,~)cb(@onVideoPan), ...
            'WindowButtonUpFcn',@(~,~)cb(@endVideoPan));
    end

    function onVideoPan()
        p = getpixelposition(axVideo,true); d = get(ancestor(axVideo,'figure'),'CurrentPoint') - vpanStart;
        r = vpanZoom; w = r(2)-r(1)+1; h = r(4)-r(3)+1;
        dpp = max(w/max(p(3),1), h/max(p(4),1));       % frame pixels per screen pixel
        x0 = round(min(max(r(1) - d(1)*dpp,1), vid.Width - w + 1));
        y0 = round(min(max(r(3) + d(2)*dpp,1), vid.Height - h + 1));
        vzoom = [x0 x0+w-1 y0 y0+h-1];
        set(axVideo,'XLim',[vzoom(1)-0.5 vzoom(2)+0.5],'YLim',[vzoom(3)-0.5 vzoom(4)+0.5]);
        if ~isempty(lastFull), setImgData(dispFrame(lastFull)); end
        drawnow limitrate;
    end

    function endVideoPan()
        vf = ancestor(axVideo,'figure');
        if isequal(vf,fig), set(fig,'WindowButtonMotionFcn',@(~,~)onHover(),'WindowButtonUpFcn','');
        else, set(vf,'WindowButtonMotionFcn','','WindowButtonUpFcn',''); end
        clearFrameCache();
    end

    % ---------------- pop-out video window ----------------------------------
    % The video axes move into their own resizable window (same zoom / pan /
    % playback); a bigger window shows the picture at higher resolution.
    function togglePopVideo()
        if ~isempty(popFig) && ishghandle(popFig), dockVideo(); else, popOutVideo(); end
    end

    function popOutVideo()
        if isempty(vid), return; end
        ss = get(0,'ScreenSize');
        w = min(1100, ss(3)-80); h = min(round(w*vid.Height/vid.Width) + 40, ss(4)-120);
        popFig = figure('Name',[APPNAME ' - video'],'NumberTitle','off','MenuBar','none', ...
            'ToolBar','none','Color',[0 0 0],'Units','pixels', ...
            'Position',[max(20,ss(3)-w-40) max(40,ss(4)-h-80) w h], ...
            'CloseRequestFcn',@(~,~)cb(@dockVideo), ...
            'SizeChangedFcn',@(~,~)cb(@onVideoResize), ...
            'WindowScrollWheelFcn',@(~,e)cb(@()onPopScroll(e)), ...
            'WindowKeyPressFcn',@(~,e)cb(@()onKey(e)));
        vidPosNorm = get(axVideo,'Position');
        set(axVideo,'Parent',popFig,'Units','normalized','Position',[0 0 1 0.93]);
        mk = @(lab,x,wd,fn,tip) uicontrol(popFig,'Style','pushbutton','String',lab, ...
            'Units','normalized','Position',[x 0.94 wd 0.05],'FontSize',11,'TooltipString',tip, ...
            'Callback',@(src,~)cb(@()barDo(src,fn)));
        mk('+',   0.01,0.05,@()videoZoomBy(1.5),  'Zoom in (or scroll on the picture)');
        mk('-',   0.07,0.05,@()videoZoomBy(1/1.5),'Zoom out');
        mk('Fit', 0.13,0.07,@resetVideoZoom,      'Whole picture (or double-click it)');
        mk('Dock',0.21,0.08,@dockVideo,           'Put the video back in the main window');
        uicontrol(popFig,'Style','text','Units','normalized','Position',[0.31 0.94 0.68 0.05], ...
            'String','Scroll = zoom   |   drag = move   |   double-click = fit   |   Space = play', ...
            'BackgroundColor',[0 0 0],'ForegroundColor',[0.85 0.85 0.85],'FontSize',10, ...
            'HorizontalAlignment','left');
        set(vidBtns(ishghandle(vidBtns)),'Visible','off');
        dockMsg = uicontrol(fig,'Style','text','Units','normalized','Position',vidPosNorm, ...
            'String',{'','','Video is in its own window.','','Click "Dock video" (or close that window)', ...
            'to bring it back here.'},'BackgroundColor',[0.15 0.15 0.15], ...
            'ForegroundColor',[1 1 1],'FontSize',11);
        set(btnPop,'String','Dock video'); uistack(btnPop,'top');
        onVideoResize();
        setStatus('Video popped out - resize that window for a bigger, sharper picture.','ok');
    end

    function dockVideo()
        if isempty(popFig) || ~ishghandle(popFig), popFig = []; return; end
        set(axVideo,'Parent',fig,'Units','normalized','Position',vidPosNorm);
        delete(popFig); popFig = [];
        deleteValid(dockMsg); dockMsg = [];
        set(vidBtns(ishghandle(vidBtns)),'Visible','on');
        set(btnPop,'String','Pop out video');
        onVideoResize();
        setStatus('Video docked back in the main window.','info');
    end

    function onVideoResize()
        % New panel size -> new resolution for the picture (and the memory copy).
        if isempty(vid) || isempty(imgVideo) || ~ishghandle(imgVideo), return; end
        updateFrameStride(); clearFrameCache();
        if ~isempty(lastFull), setImgData(dispFrame(lastFull)); end
    end

    function onPopScroll(e)
        if isempty(vid) || isEditing, return; end
        cp = get(axVideo,'CurrentPoint');
        videoZoomBy(1.15^(-e.VerticalScrollCount),cp(1,1),cp(1,2));
    end

    function onScroll(e)
        if isempty(vid) || isEditing, return; end
        if overVideo()
            cp = get(axVideo,'CurrentPoint');
            videoZoomBy(1.15^(-e.VerticalScrollCount),cp(1,1),cp(1,2));
            return;
        end
        ctrl = any(ismember({'control','command'},get(fig,'CurrentModifier')));
        c = e.VerticalScrollCount;
        if ctrl, zoomAbout(pointerTime(),1.2^c); else, scrollBy(0.15*c); end
    end


    function [a,b,tk] = snapTarget(a,b,tk)
        % A highlight that is basically one transcript piece snaps to it.
        if ~isempty(tk) || isnan(a) || isempty(events), return; end
        tr = find(strcmp({events.kind},'trans'));
        if isempty(tr), return; end
        ov = min(b,[events(tr).end]) - max(a,[events(tr).start]);
        [m,i] = max(ov);
        if m <= 0, return; end
        pl = events(tr(i)).end - events(tr(i)).start;
        if m >= 0.8*(b-a) && m >= 0.6*pl
            tk = tr(i); a = events(tk).start; b = events(tk).end;
        end
    end

    function fluentFromTarget()
        % Mark the highlighted word fluent (or pick up the fluent word already
        % there) and enter group mode.
        [a,b,tk] = currentTarget();
        if isnan(a)
            setStatus('Highlight a word in the transcript first.','warn'); return;
        end
        [a,b,tk] = snapTarget(a,b,tk);
        j = hasAnnot(a,b,'fluent','');
        if ~isempty(j), enterGroup(events(j).uid); return; end
        ov = overlapIdx(a,b,{'fluent'},[]);
        if ~isempty(ov)
            % the highlight is inside a fluent word: work on that word, and
            % treat the highlight as one of its disfluent parts
            enterGroup(events(ov(1)).uid);
            selStart = a; selEnd = b; currentEvent = [];
            updateSelectionGraphics(); redrawEvents();
            openEventPicker(); return;
        end
        pushUndo();
        k = addRaw(a,b,tk,'fluent','');
        if isempty(tk), events(k).transcript = wordsUnder(a,b); end
        relinkAll();
        eventsChanged();
        enterGroup(events(k).uid);
    end

    function w = wordsUnder(a,b)
        tr = find(strcmp({events.kind},'trans'));
        if isempty(tr), w = ''; return; end
        ov = min(b,[events(tr).end]) - max(a,[events(tr).start]);
        pl = [events(tr).end] - [events(tr).start];
        use = tr(ov > 0 & (ov > 0.3*pl | ov == max(ov)));
        [~,o] = sort([events(use).start]); use = use(o);
        w = strjoin(cellfun(@tokenWord,{events(use).transcript},'UniformOutput',false),' ');
    end

    function enterGroup(u)
        groupU = u; f = findUid(u);
        tk = findUid(events(f).tok);
        if ~isempty(tk), currentEvent = tk; else, currentEvent = f; end
        selStart = NaN; selEnd = NaN;
        updateSelectionGraphics(); redrawEvents(); syncListSelection();
        setStatus(sprintf(['FLUENT EVENT SELECTED: "%s". Now highlight each disfluent part (e.g. t-) ' ...
            'and press Enter, or W for the whole episode. Esc when done.'],events(f).transcript),'ok');
    end

    function exitGroup()
        if groupU > 0
            f = findUid(groupU); groupU = 0;
            if ~isempty(f)
                setStatus(sprintf('Finished "%s".',events(f).transcript),'ok');
            end
        end
        redrawEvents();
    end

    function wholeEpisode()
        % W: highlight the whole episode (first disfluent part to end of the word)
        if groupU <= 0
            setStatus('W works after a fluent word is selected (highlight it + Enter).','info'); return;
        end
        sp = groupSpan(groupU);
        selStart = sp(1); selEnd = sp(2); currentEvent = []; cursorTime = sp(1);
        updateSelectionGraphics(); redrawEvents();
        openEventPicker();
    end

    function openEventPicker()
        % Tick one or more disfluency / error types for the target, Enter = done.
        [a,b,tk] = currentTarget();
        if isnan(a)
            setStatus('Highlight the part first.','warn'); return;
        end
        [a,b,tk] = snapTarget(a,b,tk);
        nD = numel(DISFL); nE = numel(ERRS);
        have = false(1,nD+nE);
        for q = 1:nD, have(q) = ~isempty(hasAnnot(a,b,'disfl',DISFL(q).abbr)); end
        for q = 1:nE, have(nD+q) = ~isempty(hasAnnot(a,b,'error',ERRS(q).abbr)); end
        if ~isempty(tk), what = sprintf('"%s"',events(tk).transcript);
        else, what = sprintf('%s - %s',fmtTime(a),fmtTime(b)); end
        want = eventPicker(have,what);
        if isempty(want), setStatus('Cancelled - nothing changed.','info'); return; end
        if isequal(want,have), return; end
        applyTypes(a,b,tk,want);
    end

    function applyTypes(a,b,tk,want)
        nD = numel(DISFL);
        kinds = [repmat({'disfl'},1,nD) repmat({'error'},1,numel(ERRS))];
        cats  = [{DISFL.abbr} {ERRS.abbr}];
        keepU = 0; if ~isempty(tk), keepU = events(tk).uid; end
        pushUndo();
        del = [];
        for q = 1:numel(kinds)
            j = hasAnnot(a,b,kinds{q},cats{q});
            if ~isempty(j) && ~want(q), del(end+1) = j; end %#ok<AGROW>
        end
        events(del) = [];
        if keepU > 0, tk = findUid(keepU); end
        for q = find(want)
            if isempty(hasAnnot(a,b,kinds{q},cats{q})), addRaw(a,b,tk,kinds{q},cats{q}); end
        end
        relinkAll();
        currentEvent = []; selStart = a; selEnd = b;
        eventsChanged(); updateSelectionGraphics();
        if groupU > 0
            setStatus(sprintf(['Marked: %s. Highlight the next disfluent part + Enter, W = whole episode, ' ...
                'Esc = done.'],pieceTypes(a,b)),'ok');
        else
            setStatus(sprintf('Marked: %s (Ctrl+Z to undo).',pieceTypes(a,b)),'ok');
        end
    end

    function want = eventPicker(have,what)
        % Checkbox list of disfluency + error types. Keys R / P / B / T and
        % 1-6 tick, Enter = done, Esc = cancel. Returns [] when cancelled.
        want = [];
        nD = numel(DISFL); nE = numel(ERRS); n = nD + nE;
        rowH = 24; dw = 380; dh = 70 + (n+2)*rowH + 50;
        mp = getpixelposition(fig);
        d = figure('Name','Disfluency / error','NumberTitle','off','MenuBar','none', ...
            'ToolBar','none','Resize','off','WindowStyle','modal','Color',BG,'Units','pixels', ...
            'Position',[mp(1)+(mp(3)-dw)/2 mp(2)+(mp(4)-dh)/2 dw dh], ...
            'KeyPressFcn',@pkKey,'CloseRequestFcn',@(~,~)pkFinish(false));
        uicontrol(d,'Style','text','Units','pixels','Position',[12 dh-30 dw-24 20], ...
            'String',['What is ' what '?'],'BackgroundColor',BG,'ForegroundColor',[0 0 0], ...
            'FontWeight','bold','FontSize',11,'HorizontalAlignment','left');
        uicontrol(d,'Style','text','Units','pixels','Position',[12 dh-50 dw-24 18], ...
            'String','Tick one or more (keys R P B T and 1-6), then Enter.', ...
            'BackgroundColor',BG,'ForegroundColor',[0.25 0.25 0.25],'HorizontalAlignment','left');
        cbx = gobjects(1,n); y = dh - 78;
        uicontrol(d,'Style','text','Units','pixels','Position',[12 y dw-24 18],'String','Disfluency', ...
            'BackgroundColor',BG,'ForegroundColor',[0 0 0],'FontWeight','bold','HorizontalAlignment','left');
        for q = 1:nD
            y = y - rowH;
            cbx(q) = uicontrol(d,'Style','checkbox','Units','pixels','Position',[24 y dw-40 rowH-2], ...
                'String',sprintf('%s   %s  -  %s',upper(DISFL(q).key),DISFL(q).abbr,DISFL(q).name), ...
                'Value',have(q),'BackgroundColor',BG,'ForegroundColor',DISFL(q).color*0.75, ...
                'FontWeight','bold','KeyPressFcn',@pkKey);
        end
        y = y - rowH;
        uicontrol(d,'Style','text','Units','pixels','Position',[12 y dw-24 18],'String','Error', ...
            'BackgroundColor',BG,'ForegroundColor',[0 0 0],'FontWeight','bold','HorizontalAlignment','left');
        for q = 1:nE
            y = y - rowH;
            cbx(nD+q) = uicontrol(d,'Style','checkbox','Units','pixels','Position',[24 y dw-40 rowH-2], ...
                'String',sprintf('%d   %s  -  %s',q,ERRS(q).abbr,ERRS(q).name), ...
                'Value',have(nD+q),'BackgroundColor',BG,'ForegroundColor',ERR_COLOR*0.8, ...
                'KeyPressFcn',@pkKey);
        end
        uicontrol(d,'Style','pushbutton','Units','pixels','Position',[dw-200 12 90 30], ...
            'String','Done (Enter)','Callback',@(~,~)pkFinish(true),'KeyPressFcn',@pkKey);
        uicontrol(d,'Style','pushbutton','Units','pixels','Position',[dw-104 12 90 30], ...
            'String','Cancel (Esc)','Callback',@(~,~)pkFinish(false),'KeyPressFcn',@pkKey);
        ok = false; res = [];
        figure(d); uicontrol(cbx(1));
        uiwait(d);
        if ok, want = res; end

        function pkKey(~,ke)
            switch ke.Key
                case {'return','enter'}, pkFinish(true); return;
                case 'escape',           pkFinish(false); return;
            end
            q = find(strcmp({DISFL.key},ke.Key),1);
            if isempty(q)
                v = str2double(ke.Character);
                if isscalar(v) && ~isnan(v) && v >= 1 && v <= nE, q = nD + v; end
            end
            if ~isempty(q), set(cbx(q),'Value',1 - get(cbx(q),'Value')); end
        end
        function pkFinish(okIn)
            ok = okIn;
            if ok, res = arrayfun(@(h)get(h,'Value') == 1,cbx); end
            if ishghandle(d), delete(d); end
        end
    end

    function s = fitText(s,ax,ev,fsz)
        % Too long for its box? Show the first word + "..."; click to see it all.
        if isempty(s), return; end
        try
            px = getpixelposition(ax,true);
            w = (min(ev.end,viewEnd) - max(ev.start,viewStart))/max(viewEnd-viewStart,eps)*px(3);
        catch
            return;
        end
        n = floor(w/(0.62*fsz)) - 1;
        if numel(s) <= n, return; end
        el = char(8230);
        fw = regexp(s,'^\S+','match','once');
        if ~isempty(fw) && numel(fw) + 1 <= n, s = [fw el];
        elseif n >= 2, s = [s(1:n-1) el];
        else, s = '';
        end
    end

    function pos = dataRangeToPix(ax,x0,x1)
        % The typing box: at least 420 px wide so you can see what you write.
        axpix = getpixelposition(ax,true); xl = get(ax,'XLim');
        f0 = (x0-xl(1))/diff(xl); f1 = (x1-xl(1))/diff(xl);
        f0 = max(0,min(1,f0)); f1 = max(0,min(1,f1));
        pw = max(420,(f1-f0)*axpix(3)); pw = min(pw,axpix(3));
        px = axpix(1)+f0*axpix(3);
        px = min(px, axpix(1)+axpix(3)-pw);
        ph = max(24, min(30, axpix(4)-2));
        pos = [px axpix(2)+(axpix(4)-ph)/2 pw ph];
    end

    function editNotes()
        if isempty(currentEvent) || currentEvent > numel(events)
            setStatus(['Click an item first, then N to write its notes (or double-click a ' ...
                'highlight on the Notes strip).'],'warn'); return;
        end
        startEditing(currentEvent,'notes');
    end


    function onListSelect(~,e)
        if isempty(vid) || isempty(listMap) || isempty(e.Indices), return; end
        commitEdit();
        r = e.Indices(1,1); if r < 1 || r > numel(listMap), return; end
        if strcmp(listKind{r},'unresolved')            % open / close the unresolved section
            unresolvedOpen = ~unresolvedOpen;
            updateEventList(); return;
        end
        k = listMap(r);
        if k <= 0 || k > numel(events), return; end
        if strcmp(listKind{r},'fluent')
            % click a fluent row: open / close its disfluency / error rows
            u = events(k).uid;
            if any(expandedU == u), expandedU(expandedU == u) = [];
            else, expandedU(end+1) = u; end
            t = findUid(events(k).tok);
            if isempty(t), t = k; end
            selectEvent(t,true);
            setStatus(sprintf(['Fluent event "%s" - its disfluency / error events are listed under it. ' ...
                'Space = play, Enter = add more parts to this word.'],events(k).transcript),'info');
        else
            selectEvent(k,true);
            setStatus(sprintf('Selected %s. Enter = change its type(s), N = notes, Delete = remove.', ...
                kindLabel(k)),'info');
        end
    end

    function syncListSelection()
        % Refresh the table. Rows open / close only when you click a fluent
        % row (or the unresolved header) - selecting things never opens them.
        updateEventList();
    end

    function s = oneLine(s)
        s = strrep(strrep(char(s),sprintf('\n'),' '),sprintf('\r'),' ');
    end

    function s = fluentNotes(f)
        % notes of a fluent word = its own notes + the notes on its transcript piece
        parts = {};
        t = findUid(events(f).tok);
        if ~isempty(t) && ~isempty(events(t).notes), parts{end+1} = events(t).notes; end
        if ~isempty(events(f).notes), parts{end+1} = events(f).notes; end
        s = strjoin(parts,'; ');
    end

    function c = eventCode(k)
        % short code: F fluent, R / P / B / T disfluency, E-OS ... E-S error
        ev = events(k);
        switch ev.kind
            case 'fluent', c = 'F';
            case 'disfl'
                q = find(strcmp({DISFL.abbr},ev.cat),1);
                if isempty(q), c = upper(ev.cat); else, c = upper(DISFL(q).key); end
            case 'error',  c = ['E-' upper(regexprep(ev.cat,'^er-',''))];
            otherwise,     c = '';
        end
    end

    function C = codeTable()
        C = {'code','label','meaning','group'; 'F','fluent','fluent-accurate word (resolution)','fluent'};
        for q = 1:numel(DISFL)
            if strcmp(DISFL(q).abbr,'td'), grp = 'typical disfluency (TD)'; else, grp = 'stuttering-like disfluency (SLD)'; end
            C(end+1,:) = {upper(DISFL(q).key),DISFL(q).abbr,DISFL(q).name,grp}; %#ok<AGROW>
        end
        for q = 1:numel(ERRS)
            C(end+1,:) = {['E-' upper(regexprep(ERRS(q).abbr,'^er-',''))],ERRS(q).abbr,ERRS(q).name,'error'}; %#ok<AGROW>
        end
    end

    function C = exportCells()
        % The summary table: each fluent event, then its disfluency / error
        % events underneath; unresolved ones at the end. Times in seconds.
        hdr = {'start','end','duration','fluent_event_id','fluent_event','event_code','notes'};
        C = hdr;
        if isempty(events), return; end
        fl = find(strcmp({events.kind},'fluent'));
        [~,o] = sort([events(fl).start]); fl = fl(o);
        r3 = @(x) round(x*1000)/1000;
        for i = 1:numel(fl)
            f = fl(i);
            C(end+1,:) = {r3(events(f).start), r3(events(f).end), r3(events(f).end-events(f).start), ...
                i, events(f).transcript, 'F', fluentNotes(f)}; %#ok<AGROW>
            kids = linkedAll(events(f).uid);
            [~,o] = sort([events(kids).start]); kids = kids(o);
            for j = kids
                C(end+1,:) = {r3(events(j).start), r3(events(j).end), r3(events(j).end-events(j).start), ...
                    i, events(f).transcript, eventCode(j), events(j).notes}; %#ok<AGROW>
            end
        end
        orph = find(ismember({events.kind},{'disfl','error'}) & ...
            arrayfun(@(e) e.link <= 0 || isempty(findUid(e.link)),events));
        [~,o] = sort([events(orph).start]); orph = orph(o);
        for j = orph
            C(end+1,:) = {r3(events(j).start), r3(events(j).end), r3(events(j).end-events(j).start), ...
                '', '(unresolved)', eventCode(j), events(j).notes}; %#ok<AGROW>
        end
    end

    function C = readSavedTable(f)
        % Reopen: the annotator_data sheet (new files) or the old layouts.
        [d,b,e] = fileparts(f);
        if strcmpi(e,'.csv')
            b = regexprep(b,'_(annotator|codes)$','');
            alt = fullfile(d,[b '_annotator.csv']);
            if exist(alt,'file'), C = readAnnotCells(alt,'events'); return; end
            C = readAnnotCells(f,'events'); return;
        end
        try
            C = readAnnotCells(f,'annotator_data'); return;
        catch
        end
        C = readAnnotCells(f,'events');
    end

    % ======================== unchanged helpers ==============================

    function cb(fn)
        try
            fn();
        catch err
            loc = '';
            frm = [];
            for si = 1:numel(err.stack)
                if strcmp(err.stack(si).name,mfilename) || ...
                        startsWith(err.stack(si).name,[mfilename '/'])
                    frm = err.stack(si); break;
                end
            end
            if isempty(frm) && ~isempty(err.stack), frm = err.stack(1); end
            if ~isempty(frm)
                loc = sprintf('\n\n(in %s, line %d)',frm.name,frm.line);
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
        % Toggling Enable removes keyboard focus from the button, so keys go
        % to the figure and do not re-click it.
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
        text(sax,0.5,0.66,['Fluency, disfluency and error annotation  |  v' APPVER], ...
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

    function uiLoadVideo()
        if ~confirmUnsaved('loading another video'), return; end
        [fn,fp] = uigetfile( ...
            {'*.mp4;*.mov;*.avi;*.m4v;*.mkv;*.wmv','Video files';'*.*','All files'}, ...
            'Select a video file');
        if isequal(fn,0), return; end
        loadVideo(fullfile(fp,fn));
    end

    function loadVideo(path)
        path = originalOf(path);             % picked a _smooth copy? use its original
        vpath = smoothCopy(path);            % fast-decoding copy for the picture
        wb = waitbar(0,'Opening video file...','Name','Loading video');
        wbClean = onCleanup(@()deleteValid(wb)); %#ok<NASGU>
        try
            v = VideoReader(vpath);
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
        events = EV0; nextUid = 1;
        currentEvent = []; selStart = NaN; selEnd = NaN; cursorTime = 0;
        pendingNew = []; groupU = 0;
        vzoom = []; lastFull = []; clearFrameCache();
        undoStack = {}; dirty = false;
        viewStart = 0; viewEnd = min(dur,10);
        buildPlayer();
        try
            % tick about twice per video frame so no frame is skipped or doubled
            set(playTimer,'Period',max(0.010,round(1000*0.5/max(vid.FrameRate,1))/1000));
        catch
        end

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
        updateFrameStride();
        imgVideo = image(axVideo,'XData',[1 max(vid.Width,2)],'YData',[1 max(vid.Height,2)], ...
            'CData',dispFrame(frame0));
        set(axVideo,'YDir','reverse'); axis(axVideo,'image','off');
        set(imgVideo,'PickableParts','none','HitTest','off');
        lastFull = frame0;

        set(videoDependent,'Enable','on');
        updateEventList(); updateTitle();
        waitbar(1,wb,'Done');
        refreshView(); showFrameAt(0);
        [~,nm,ext] = fileparts(path);
        if strcmp(vpath,path), how = ''; else, how = ', smooth copy'; end
        setStatus(sprintf(['Loaded %s%s  (%s, %d Hz audio, %.0f fps%s).  Drag to select speech, ' ...
            'then double-click it (or Enter) to type the transcript.'], ...
            nm,ext,fmtTime(dur),fs,vid.FrameRate,how),'ok');
    end

    function p = originalOf(p)
        % If a _smooth copy was picked, find the original next to it (it has the audio).
        [d,b] = fileparts(p);
        if ~endsWith(b,SMOOTH_SUFFIX), return; end
        b0 = b(1:end-numel(SMOOTH_SUFFIX));
        for ex = {'.MOV','.mov','.mp4','.MP4','.m4v','.avi','.mkv','.wmv'}
            c = fullfile(d,[b0 ex{1}]);
            if exist(c,'file'), p = c; return; end
        end
    end

    function vpath = smoothCopy(path)
        % Make (once) and return a copy that MATLAB decodes fast: H.264,
        % at most SMOOTH_HEIGHT tall, constant frame rate, no audio.
        % Falls back to the original file whenever anything goes wrong.
        vpath = path;
        if ~SMOOTH_ENABLE, return; end
        [d,b,e] = fileparts(path);
        out = fullfile(d,[b SMOOTH_SUFFIX '.mp4']);
        if exist(out,'file')
            try
                vo = VideoReader(out);
                if vo.Duration > 0, vpath = out; return; end
            catch
            end
            try, delete(out); catch, end            % broken copy: make it again
        end
        fps = 30; need = ~any(strcmpi(e,{'.mp4','.m4v'}));
        try
            v0 = VideoReader(path); fps = v0.FrameRate;
            need = need || v0.Height > SMOOTH_HEIGHT;
        catch
            need = true;
        end
        if ~need, return; end
        ff = findFfmpeg();
        if isempty(ff)
            msgbox(sprintf(['For smoother video, install ffmpeg once and reload the video:\n\n' ...
                'Mac (Terminal):  brew install ffmpeg\n\nUsing the original file for now.']), ...
                'Smoother playback','help','replace');
            return;
        end
        auto = false;
        try, auto = ispref('DisfluencyAnnotator','autoSmooth') && ...
                getpref('DisfluencyAnnotator','autoSmooth'); catch, end
        if ~auto
            c = questdlg(sprintf(['This video is large or in a format MATLAB decodes slowly.\n\n' ...
                'Make a smooth playback copy (%dp, constant frame rate)? It is made once ' ...
                '(about as long as the video, often less) and saved next to it as\n%s\n\n' ...
                'Audio and annotation times still come from the original.'], ...
                SMOOTH_HEIGHT,[b SMOOTH_SUFFIX '.mp4']), ...
                'Smoother playback','Make copy','Always','Use original','Make copy');
            switch c
                case 'Always', try, setpref('DisfluencyAnnotator','autoSmooth',true); catch, end
                case 'Make copy'
                otherwise, return;
            end
        end
        fpsOut = min(60,max(1,round(fps)));
        cmd = sprintf(['"%s" -y -loglevel error -i "%s" -an -vf "scale=-2:min(%d\\,ih)" ' ...
            '-r %d -c:v libx264 -preset veryfast -crf 20 -pix_fmt yuv420p ' ...
            '-movflags +faststart "%s"'],ff,path,SMOOTH_HEIGHT,fpsOut,out);
        wb2 = waitbar(0.5,sprintf('Making smooth playback copy (one time only)...'), ...
            'Name','Smooth copy');
        [st,msg] = system(cmd);
        deleteValid(wb2);
        if st == 0 && exist(out,'file')
            vpath = out;
        else
            try, if exist(out,'file'), delete(out); end; catch, end
            warndlg(sprintf('Could not make the smooth copy - using the original.\n\n%s', ...
                strtrim(msg)),'Smooth copy');
        end
    end

    function ff = findFfmpeg()
        % MATLAB's PATH often misses Homebrew, so check the usual places first.
        ff = '';
        for c = {'/opt/homebrew/bin/ffmpeg','/usr/local/bin/ffmpeg','/usr/bin/ffmpeg'}
            if exist(c{1},'file'), ff = c{1}; return; end
        end
        try
            if ispc, [s,r] = system('where ffmpeg'); else, [s,r] = system('command -v ffmpeg'); end
            if s == 0
                r = strtrim(strsplit(strtrim(r),newline));
                if ~isempty(r{1}) && exist(r{1},'file'), ff = r{1}; end
            end
        catch
        end
    end

    function refreshView()
        if isempty(vid), return; end
        viewStart = max(0,viewStart); viewEnd = min(dur,viewEnd);
        if viewEnd - viewStart < MINVIEW, viewEnd = min(dur,viewStart+MINVIEW); end
        set([axSpec axWave axTrans axDE axNotes],'XLim',[viewStart viewEnd]);
        set([axSpec axWave axTrans axDE],'XTick',[]); set(axNotes,'XTickMode','auto');
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
        if numel(seg) < winlen, return; end
        Ncol = 700;
        hop  = max(1, floor((numel(seg)-winlen)/Ncol));
        nfft = 2^nextpow2(max(winlen,512));
        w    = 0.5 - 0.5*cos(2*pi*(0:winlen-1)'/(winlen-1));   % Hann window (no toolbox)
        st   = 1:hop:(numel(seg)-winlen+1);
        if isempty(st), return; end
        nf = floor(nfft/2)+1;
        % all frames in one FFT call (much faster than a loop)
        idx = st + (0:winlen-1)';
        X = fft(seg(idx).*w, nfft);
        P = 20*log10(abs(X(1:nf,:))+eps);
        F = (0:nf-1)'*(fs/nfft);
        Tc = (st-1+winlen/2)/fs + viewStart;
        fmask = F <= maxFreq; if ~any(fmask), fmask = true(size(F)); end
        hSpecImg = imagesc(axSpec, Tc, F(fmask), P(fmask,:));
        set(axSpec,'YDir','normal','YLim',[0 maxFreq],'XLim',[viewStart viewEnd]);
        mx = max(max(P(fmask,:)));
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

    function h = drawHandles(ax,t0,t1)
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

    function onOverviewDown()
        if isempty(vid), return; end
        commitEdit();
        cp = get(axOverview,'CurrentPoint'); t = clampT(cp(1,1));
        ovMode = ovHit(t);
        if isempty(ovMode)
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
        set(box,'UserData',get(box,'String'),'String','','Enable','on'); uicontrol(box);
        setStatus(['Type a time - e.g. 14 = 14 s, 1:05 = 1 min 5 s, 2m30 also works - ' ...
            'then press Enter or click Go (Esc cancels).'],'info');
    end

    function syncTimeFields(a,b)
        if isempty(edtFrom) || ~ishghandle(edtFrom), return; end
        if ~isequal(typingBox,edtFrom), set(edtFrom,'String',fmtTime(a)); end
        if ~isequal(typingBox,edtTo),   set(edtTo,'String',fmtTime(b));   end
    end

    function applyTimeFields(changed)
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

    % ---------------- zoom / scroll ------------------------------------------

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

    function fullView()
        setView(0,dur);
    end

    function zoomToSelection()
        if ~hasSel()
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

    % ---------------- mouse on the timelines --------------------------------
    % Press on a black edge handle -> resize.  Press inside an event -> drag to
    % move it (a plain click selects it).  Press elsewhere (or Shift + press)
    % -> drag a new selection.

    function revertLast(msg)
        if isempty(undoStack), return; end
        events = undoStack{end}; undoStack(end) = [];
        if ~isempty(currentEvent) && currentEvent <= numel(events)
            selStart = events(currentEvent).start; selEnd = events(currentEvent).end;
        end
        eventsChanged(); updateSelectionGraphics();
        setStatus(msg,'warn'); beep;
    end

    function edge = hitEdge(ax,t)
        edge = '';
        if ~hasSel(), return; end
        tol = EDGE_PX*timePerPixel(ax);
        dL = abs(t-selStart); dR = abs(t-selEnd);
        if min(dL,dR) > tol, return; end
        if dL <= dR, edge = 'left'; else, edge = 'right'; end
    end

    function tpp = timePerPixel(ax)
        p = getpixelposition(ax,true); tpp = (viewEnd-viewStart)/max(p(3),1);
    end

    function tf = justCommitted()
        % A text box handles Enter / Esc itself; the figure may see the same
        % key press a moment later. Ignore it.
        tf = (now*86400 - lastKeyCommit) < 0.4;
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

    % ---------------- adding events -----------------------------------------

    function q = errorPicker(titleStr)
        % Clickable list of error types; 1-6 also pick, Esc / Cancel cancels.
        q = [];
        n = numel(ERRS); bh = 30; gap = 6; dw = 330;
        dh = 34 + (n+1)*(bh+gap) + gap;
        mp = getpixelposition(fig);
        d = figure('Name',titleStr,'NumberTitle','off','MenuBar','none','ToolBar','none', ...
            'Resize','off','WindowStyle','modal','Color',BG,'Units','pixels', ...
            'Position',[mp(1)+(mp(3)-dw)/2 mp(2)+(mp(4)-dh)/2 dw dh], ...
            'KeyPressFcn',@dKey,'CloseRequestFcn',@(~,~)finish(0));
        uicontrol(d,'Style','text','String','Choose the error type  (keys 1-6, Esc = cancel)', ...
            'Units','pixels','Position',[10 dh-28 dw-20 20],'BackgroundColor',BG, ...
            'ForegroundColor',[0 0 0],'FontWeight','bold');
        for ii = 1:n
            yy = dh - 34 - ii*(bh+gap) + gap;
            uicontrol(d,'Style','pushbutton','Units','pixels','Position',[12 yy dw-24 bh], ...
                'String',sprintf('%d.   %s   -   %s',ii,ERRS(ii).abbr,ERRS(ii).name), ...
                'Callback',@(~,~)finish(ii),'KeyPressFcn',@dKey);
        end
        yy = dh - 34 - (n+1)*(bh+gap) + gap;
        uicontrol(d,'Style','pushbutton','Units','pixels','Position',[12 yy dw-24 bh], ...
            'String','Cancel  (Esc)','Callback',@(~,~)finish(0),'KeyPressFcn',@dKey);
        choice = 0;
        uiwait(d);
        if choice > 0, q = choice; end

        function dKey(~,ke)
            if strcmp(ke.Key,'escape'), finish(0); return; end
            v = str2double(ke.Character);
            if isscalar(v) && ~isnan(v) && v >= 1 && v <= n, finish(v); end
        end
        function finish(v)
            choice = v;
            if ishghandle(d), delete(d); end
        end
    end

    function eventsChanged()
        dirty = true; updateTitle(); updateEventList(); redrawEvents();
    end

    function k = findUid(u)
        k = [];
        if isempty(events) || isempty(u) || u <= 0, return; end
        k = find([events.uid] == u,1);
    end

    function clearGlobal(k)
        % Global start/stop no longer match once an event is moved / resized.
        events(k).gStart = NaN; events(k).gStop = NaN;
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

    function s = eventLabel(k)
        s = sprintf('%s (%s - %s)',kindLabel(k),fmtTime(events(k).start),fmtTime(events(k).end));
    end

    function s = disflName(abbr)
        q = find(strcmp({DISFL.abbr},abbr),1);
        if isempty(q), s = 'unknown type'; else, s = DISFL(q).name; end
    end

    function s = errName(abbr)
        q = find(strcmp({ERRS.abbr},abbr),1);
        if isempty(q), s = 'unknown type'; else, s = ERRS(q).name; end
    end

    % ---------------- annotation text editing -------------------------------

    function commitEdit()
        % Let a pending edit-box callback run first (macOS commits the text
        % on focus loss), then save whatever is still being typed.
        drawnow;
        if isEditing && ~isempty(typingBox) && ishghandle(typingBox)
            commitBox(typingBox);
        end
        isEditing = false; typingBox = [];
    end

    function pushUndo()
        undoStack{end+1} = events;
        if numel(undoStack) > MAXUNDO, undoStack(1) = []; end
    end

    function buildPlayer()
        try, if ~isempty(mainPlayer), stop(mainPlayer); end; catch, end
        mainPlayer = [];
        try
            if audiodevinfo(0) < 1
                errordlg(['MATLAB cannot find an audio output device, so playback ' ...
                    'will be silent. Check your sound output and restart MATLAB.'], ...
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

    function restartLoop()
        % Jump back to the start of the bounded range and keep playing.
        s0 = max(1,round(playLoopA*fs)+1); s1 = min(numel(audio),round(playEndTime*fs));
        if s1 <= s0, stopPlay(); return; end
        try
            stop(mainPlayer);
            try
                vid.CurrentTime = max(0, min(playLoopA, vid.Duration - 1/max(vid.FrameRate,1)));
            catch
            end
            play(mainPlayer,[s0 s1]);
        catch
            stopPlay(); return;
        end
        loopPass = loopPass + 1;
        if ~isempty(playGfx) && all(ishghandle(playGfx))
            set(playGfx(1),'XData',[playLoopA playLoopA]); set(playGfx(2),'XData',[playLoopA playLoopA]);
        end
        setStatus(sprintf('Looping %s - %s (pass %d).   Space = pause, L = stop looping.', ...
            fmtTime(playLoopA),fmtTime(playEndTime),loopPass),'info');
    end

    function t = currentPlayTime(player)
        t = playStartTime + double(player.CurrentSample - 1)/fs;
    end

    function pausePlay()
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

    function v = cellOr(C,i,c)
        if c == 0 || c > size(C,2), v = ''; else, v = C{i,c}; end
    end

    function x = firstOr0(x)
        if isempty(x), x = 0; end
    end

    function s = rowList(v)
        s = strjoin(arrayfun(@(x)sprintf('%d',x),v(1:min(6,numel(v))),'UniformOutput',false),', ');
        if numel(v) > 6, s = [s sprintf(', ... %d total',numel(v))]; end
    end

    function v = numOrEmpty(x)
        if isempty(x) || isnan(x), v = ''; else, v = x; end
    end

    function s = orNone(s,alt)
        if isempty(s), s = alt; end
    end

    function tf = hasSel()
        tf = ~isnan(selStart) && ~isnan(selEnd) && selEnd > selStart;
    end

    function C = readAnnotCells(f,sheet)
        % Read an annotation table as a cell array (row 1 = headers).
        [~,~,ext] = fileparts(f); ext = lower(ext); isEvents = strcmp(sheet,'events');
        msgs = {};
        try
            if isEvents, C = readcell(f); else, C = readcell(f,'Sheet',sheet); end
            return;
        catch err
            msgs{end+1} = ['readcell: ' err.message];
        end
        try
            if isEvents, T = readtable(f,'TextType','char');
            else,        T = readtable(f,'Sheet',sheet,'TextType','char'); end
            C = [T.Properties.VariableNames; table2cell(T)]; return;
        catch err
            msgs{end+1} = ['readtable: ' err.message];
        end
        try
            switch ext
                case '.csv'
                    if ~isEvents, error('CSV files have no time-frame sheet.'); end
                    C = readCsvManual(f);
                case '.xlsx'
                    C = readXlsxManual(f,sheet);
                otherwise
                    error('No built-in reader for %s files - save it as .xlsx or .csv.',ext);
            end
            return;
        catch err
            msgs{end+1} = ['built-in reader: ' err.message];
        end
        error('%s',strjoin(regexprep(msgs,'<[^>]*>',''),sprintf('\n')));
    end

    function C = readCsvManual(f)
        txt = fileread(f);
        if numel(txt) >= 1 && double(txt(1)) == 65279, txt = txt(2:end); end   % BOM
        rows = {}; field = ''; row = {}; inQ = false; k = 1; n = numel(txt);
        while k <= n
            ch = txt(k);
            if inQ
                if ch == '"'
                    if k < n && txt(k+1) == '"', field(end+1) = '"'; k = k + 1; %#ok<AGROW>
                    else, inQ = false; end
                else
                    field(end+1) = ch; %#ok<AGROW>
                end
            elseif ch == '"'
                inQ = true;
            elseif ch == ','
                row{end+1} = field; field = ''; %#ok<AGROW>
            elseif ch == 10 || ch == 13
                if ch == 13 && k < n && txt(k+1) == 10, k = k + 1; end
                row{end+1} = field; field = ''; %#ok<AGROW>
                if ~(numel(row) == 1 && isempty(row{1})), rows{end+1} = row; end %#ok<AGROW>
                row = {};
            else
                field(end+1) = ch; %#ok<AGROW>
            end
            k = k + 1;
        end
        if ~isempty(field) || ~isempty(row), row{end+1} = field; rows{end+1} = row; end
        if isempty(rows), C = cell(0,0); return; end
        nc = max(cellfun(@numel,rows)); C = repmat({''},numel(rows),nc);
        for r = 1:numel(rows), C(r,1:numel(rows{r})) = rows{r}; end
        for r = 2:size(C,1)
            for q = 1:nc
                v = str2double(C{r,q}); if ~isnan(v), C{r,q} = v; end
            end
        end
    end

    function C = readXlsxManual(f,sheet)
        % Minimal .xlsx reader: unzip and parse the sheet XML with REGEXP.
        tmp = tempname; unzip(f,tmp);
        cleanTmp = onCleanup(@()rmdir(tmp,'s')); %#ok<NASGU>
        wsDir = fullfile(tmp,'xl','worksheets');
        if strcmp(sheet,'events'), sheetFile = fullfile(wsDir,'sheet1.xml');
        else
            wb = fileread(fullfile(tmp,'xl','workbook.xml'));
            names = regexp(wb,'<sheet[^>]*name="([^"]*)"','tokens');
            names = cellfun(@(c)c{1},names,'UniformOutput',false);
            idx = find(strcmpi(names,sheet),1);
            if isempty(idx), error('Sheet "%s" not found.',sheet); end
            sheetFile = fullfile(wsDir,sprintf('sheet%d.xml',idx));
        end
        if ~exist(sheetFile,'file'), error('Worksheet not found inside the .xlsx file.'); end
        shared = {};
        ssFile = fullfile(tmp,'xl','sharedStrings.xml');
        if exist(ssFile,'file')
            ss = fileread(ssFile);
            si = regexp(ss,'<si>(.*?)</si>','tokens');
            shared = cell(1,numel(si));
            for q = 1:numel(si)
                tt = regexp(si{q}{1},'<t[^>]*>(.*?)</t>','tokens');
                shared{q} = xmlUnescape(strjoin(cellfun(@(c)c{1},tt,'UniformOutput',false),''));
            end
        end
        x = fileread(sheetFile);
        cells = regexp(x,'<c r="([A-Z]+)(\d+)"([^>]*?)(/>|>(.*?)</c>)','tokens');
        C = {};
        for q = 1:numel(cells)
            colL = cells{q}{1}; r = str2double(cells{q}{2}); attrs = cells{q}{3};
            body = ''; if numel(cells{q}) >= 5, body = cells{q}{5}; end
            cIdx = 0; for L = colL, cIdx = cIdx*26 + (double(L)-64); end
            v = regexp(body,'<v>(.*?)</v>','tokens','once');
            tAttr = regexp(attrs,'t="([^"]*)"','tokens','once');
            if ~isempty(tAttr) && strcmp(tAttr{1},'s') && ~isempty(v)
                val = shared{str2double(v{1})+1};
            elseif ~isempty(tAttr) && strcmp(tAttr{1},'inlineStr')
                tt = regexp(body,'<t[^>]*>(.*?)</t>','tokens');
                val = xmlUnescape(strjoin(cellfun(@(c)c{1},tt,'UniformOutput',false),''));
            elseif ~isempty(tAttr) && strcmp(tAttr{1},'str') && ~isempty(v)
                val = xmlUnescape(v{1});
            elseif ~isempty(v)
                val = str2double(v{1});
            else
                val = '';
            end
            C{r,cIdx} = val; %#ok<AGROW>
        end
        C(cellfun(@isempty,C)) = {''};
    end

    function s = xmlUnescape(s)
        s = strrep(s,'&lt;','<'); s = strrep(s,'&gt;','>'); s = strrep(s,'&quot;','"');
        s = strrep(s,'&apos;',''''); s = strrep(s,'&amp;','&');
    end

    function writeCsvManual(out,header,data)
        fid = fopen(out,'w','n','UTF-8');
        if fid < 0, error('Cannot open %s for writing.',out); end
        closer = onCleanup(@()fclose(fid)); %#ok<NASGU>
        fprintf(fid,'%s\n',strjoin(header,','));
        for r = 1:size(data,1)
            parts = cell(1,size(data,2));
            for q = 1:size(data,2)
                v = data{r,q};
                if isnumeric(v)
                    if isempty(v) || isnan(v), parts{q} = ''; else, parts{q} = sprintf('%.10g',v); end
                else
                    parts{q} = ['"' strrep(char(v),'"','""') '"'];
                end
            end
            fprintf(fid,'%s\n',strjoin(parts,','));
        end
    end

    function [isDup,dupUid] = compareToCurrent(newEv)
        isDup = false(1,numel(newEv)); dupUid = zeros(1,numel(newEv));
        es = [events.start]; ee = [events.end];
        key = strcat({events.kind},'|',{events.cat});
        for i = 1:numel(newEv)
            d = find(abs(es-newEv(i).start) < 1e-3 & abs(ee-newEv(i).end) < 1e-3 & ...
                strcmp(key,[newEv(i).kind '|' newEv(i).cat]),1);
            if ~isempty(d), isDup(i) = true; dupUid(i) = events(d).uid; end
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

    function deleteValid(h)
        if isempty(h), return; end
        h = h(ishghandle(h)); if ~isempty(h), delete(h); end
    end

    function stackEach(objs,where)
        objs = objs(ishghandle(objs));
        for o = reshape(objs,1,[])
            uistack(o,where);
        end
    end

    function x = toNum(v)
        if iscell(v), if isempty(v), x = NaN; return; end; v = v{1}; end
        if isempty(v) || isa(v,'missing'), x = NaN; return; end
        if isnumeric(v) || islogical(v), x = double(v(1)); return; end
        try, x = str2double(string(v)); catch, x = NaN; end
        if isempty(x) || ismissing(x), x = NaN; end
    end

    function s = toStr(v)
        if iscell(v), if isempty(v), s = ''; return; end; v = v{1}; end
        if isempty(v) || isa(v,'missing'), s = ''; return; end
        if ischar(v), s = strtrim(v); return; end
        if isnumeric(v) || islogical(v)
            if isnan(double(v(1))), s = ''; else, s = num2str(v(1)); end
            return;
        end
        try
            if any(ismissing(v)), s = ''; return; end
            s = strtrim(char(string(v)));
        catch
            s = '';
        end
        if strcmpi(s,'NaN') || strcmp(s,'<missing>'), s = ''; end
    end

    function onClose()
        try
            if ~confirmUnsaved('closing'), return; end
        catch
        end
        try, stop(playTimer); end %#ok<TRYNC>
        try, delete(playTimer); end %#ok<TRYNC>
        try, if ~isempty(mainPlayer), stop(mainPlayer); end; catch, end
        try, if ~isempty(popFig) && ishghandle(popFig), delete(popFig); end; catch, end
        delete(fig);
    end
end