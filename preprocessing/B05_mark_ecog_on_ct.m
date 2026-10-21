function mark_ecog_on_ct()
% mark_ecog_on_ct
% Interactive GUI for manually marking ECoG electrode locations on a CT scan (.nii) in 3-D.
%
% OVERVIEW:
%   Displays a NIfTI CT volume in three orthogonal views (Sagittal/X, Coronal/Y, Axial/Z)
%   and lets the user click to place electrode markers. Electrode positions are stored in
%   both voxel-space and millimetre coordinates and can be saved/loaded as .tsv files.
%
% ELECTRODE NAMING CONVENTION:
%   Strip 1 (posterior): labels 101-163  (3 columns x 21 rows; base 100)
%   Strip 2 (anterior):  labels 201-263  (3 columns x 21 rows; base 200)
%   Formula: label = strip_base + (col_idx-1)*21 + row_idx
%
% CONTROLS (also listed in the GUI):
%   Mouse wheel scroll         - shift slice plane in the hovered view (1 slice/tick)
%   Shift + mouse wheel scroll - shift slice plane in the hovered view (4 slices/tick)
%   e   - activate mark mode for the next unmarked electrode
%   d   - toggle delete-marker mode
%   o   - toggle out-of-plane electrode display
%   c   - toggle crosshairs
%   m   - maximise / restore the hovered view
%   Ctrl+O - load electrode locations from TSV
%   Ctrl+S - save electrode locations to TSV

    %% ======================================================================
    %% 1.  SUBJECT & FILE SETUP
    %% ======================================================================

    options.subject_id = 'DM10';    % default subject identifier

    % Prompt user to confirm or update the subject ID
    user_answer = inputdlg({'Enter Subject ID:'}, 'Subject Setup', [1 35], {options.subject_id});
    if isempty(user_answer)
        disp('Cancelled by user.');
        return;
    end
    options.subject_id = user_answer{1};

    % Attempt to auto-detect the root data directory; fall back to a picker dialog
    if     exist('Y:\',             'dir'),  options.root_dir = 'Y:\';
    elseif exist('/Volumes/Nexus4', 'dir'),  options.root_dir = '/Volumes/Nexus4';
    else
        options.root_dir = uigetdir('', 'Select root directory (should contain folder ''DBS'')');
        if options.root_dir == 0
            disp('No root directory selected. Exiting.');
            return;
        end
    end

    %% ======================================================================
    %% 2.  SHARED APPLICATION STATE
    %% ======================================================================
    % These variables are visible to — and modifiable by — all nested functions
    % below, via MATLAB's nested-function workspace-sharing mechanism.

    brain_volume      = [];         % 3-D double array of CT voxel intensities  (X x Y x Z)
    nifti_info        = [];         % struct from niftiinfo(); holds the affine transform
    volume_size       = [1 1 1];    % [nX, nY, nZ]  dimensions of the loaded volume
    current_slices    = [1 1 1];    % [iX, iY, iZ]  indices of the currently displayed planes

    active_electrode  = NaN;        % numeric label of electrode being placed (NaN = none active)
    delete_mode       = false;      % when true, clicking near an electrode removes it
    show_out_of_plane = true;       % when true, electrodes outside the current slice appear faintly
    show_crosshairs   = true;       % when true, intersection crosshair lines are drawn

    hovered_view_name = '';         % full name ('Sagittal (X)' etc.) of the view under the cursor
    maximized_view    = '';         % name of the currently maximised view  ('' = all equal)
    hovered_electrode = NaN;        % in delete mode: label of electrode nearest the cursor

    unsaved_electrode_changes = false;  % true whenever electrode data changes without saving

    % Handles for the temporary hover-ring graphics drawn in delete mode.
    % Stored per-view so they can be removed on the next hover update.
    hover_highlight_handles = struct('sag', [], 'cor', [], 'axi', []);

    % Per-view display orientation flags.  These are DISPLAY-ONLY — they never alter voxel
    % coordinates or click-to-voxel mapping.  MATLAB's CurrentPoint always returns values in
    % the axes data-coordinate space regardless of XDir / YDir settings, so electrode
    % placement and scatter overlays remain correct regardless of which flips are active.
    %
    %   _ud flags  control the vertical   display direction (YDir):
    %              false -> YDir = 'reverse' (row 1 at top;    imagesc default)
    %              true  -> YDir = 'normal'  (row 1 at bottom; image flipped vertically)
    %
    %   _lr flags  control the horizontal display direction (XDir):
    %              false -> XDir = 'normal'  (col 1 at left;   imagesc default)
    %              true  -> XDir = 'reverse' (col 1 at right;  image flipped horizontally)
    flip_sagittal_view_ud = false;
    flip_coronal_view_ud  = false;
    flip_axial_view_ud    = false;
    flip_sagittal_view_lr = false;
    flip_coronal_view_lr  = false;
    flip_axial_view_lr    = false;

    %% ======================================================================
    %% 3.  ELECTRODE TABLE INITIALISATION
    %% ======================================================================

    electrode_labels = [];
    for strip_idx = 1:2
        strip_base = strip_idx * 100;           % 100 for strip 1 (posterior), 200 for strip 2 (anterior)
        for col_idx = 1:3
            for row_idx = 1:21
                electrode_labels = [electrode_labels; ...
                    strip_base + (col_idx - 1) * 21 + row_idx]; %#ok<AGROW>
            end
        end
    end

    num_electrodes   = numel(electrode_labels);
    electrode_colors = jet(num_electrodes);     % one unique jet colour per electrode

    % Master data table — one row per electrode
    electrodes = table( ...
        electrode_labels, ...
        nan(num_electrodes,1), nan(num_electrodes,1), nan(num_electrodes,1), ...  % voxel_x/y/z
        nan(num_electrodes,1), nan(num_electrodes,1), nan(num_electrodes,1), ...  % mm_x/y/z
        electrode_colors, ...
        'VariableNames', {'name','voxel_x','voxel_y','voxel_z','mm_x','mm_y','mm_z','color'});

    %% ======================================================================
    %% 4.  MAIN FIGURE & TOP-LEVEL LAYOUT
    %% ======================================================================

    main_fig = uifigure('Name', 'Mark ECoG on CT', 'Position', [100 100 1500 950]);
    main_fig.WindowKeyPressFcn     = @handle_keypress;
    main_fig.WindowScrollWheelFcn  = @handle_scroll;
    main_fig.WindowButtonMotionFcn = @track_mouse_hover;
    main_fig.WindowButtonDownFcn   = @handle_fig_click;

    % Two-column layout: left = electrode list panel, right = brain viewer panel
    top_layout = uigridlayout(main_fig, [1 2]);
    top_layout.ColumnWidth = {'1.2x', '2.5x'};

    %% ======================================================================
    %% 5.  LEFT PANEL — ELECTRODE STRIPS & CONTROLS
    %% ======================================================================

    electrode_panel        = uipanel(top_layout);
    electrode_panel_layout = uigridlayout(electrode_panel, [3 1]);
    electrode_panel_layout.RowHeight = {'fit', '1x', 'fit'};
    electrode_panel_layout.Padding   = [5 2 5 2];

    % ---- Top row: file-operation buttons ----
    file_btn_row = uigridlayout(electrode_panel_layout, [1 3]);
    file_btn_row.ColumnWidth = {'1x', '1x', '1x'};
    file_btn_row.Padding     = [0 0 0 0];
    uibutton(file_btn_row, 'Text', 'Load CT', ...
        'ButtonPushedFcn', @(~,~) load_brain_volume());
    uibutton(file_btn_row, 'Text', 'Load Electrodes (Ctrl+O)', ...
        'ButtonPushedFcn', @(~,~) load_electrode_locations());
    uibutton(file_btn_row, 'Text', 'Save Electrodes (Ctrl+S)', ...
        'ButtonPushedFcn', @(~,~) save_electrode_locations());

    % ---- Middle: scrollable electrode grid ----
    % 7 display columns:
    %   cols 1-3  = Strip 2 (anterior)   [col 1 = outermost, col 3 = innermost]
    %   col  4    = dark visual divider
    %   cols 5-7  = Strip 1 (posterior)  [col 5 = outermost, col 7 = innermost]
    strip_scroll_panel = uipanel(electrode_panel_layout, 'Scrollable', 'on');
    strip_grid         = uigridlayout(strip_scroll_panel, [22 7]);
    strip_grid.ColumnWidth = {'1x','1x','1x', 5, '1x','1x','1x'};
    strip_grid.RowHeight   = [{18}, repmat({24}, 1, 21)];   % row 1 = header; rows 2-22 = electrodes

    % Strip header labels (row 1 of the grid)
    lbl_strip2 = uilabel(strip_grid, 'Text', 'Strip 2 (anterior)', ...
        'HorizontalAlignment','center', 'FontWeight','bold');
    lbl_strip2.Layout.Row    = 1;
    lbl_strip2.Layout.Column = [1 3];

    lbl_strip1 = uilabel(strip_grid, 'Text', 'Strip 1 (posterior)', ...
        'HorizontalAlignment','center', 'FontWeight','bold');
    lbl_strip1.Layout.Row    = 1;
    lbl_strip1.Layout.Column = [5 7];

    % Dark divider bar between the two strips (column 4, all data rows)
    for divider_row = 2:22
        div_bar = uipanel(strip_grid, 'BackgroundColor',[0.2 0.2 0.2], 'BorderType','none');
        div_bar.Layout.Row    = divider_row;
        div_bar.Layout.Column = 4;
    end

    % Map display (grid_row x display_column) -> electrode label.
    % Strip 1 (posterior) uses display columns 5-7:
    %   col 7 = column A (labels 101-121), col 6 = col B (122-142), col 5 = col C (143-163)
    % Strip 2 (anterior) uses display columns 1-3:
    %   col 3 = column A (labels 201-221), col 2 = col B (222-242), col 1 = col C (243-263)
    electrode_grid_map = nan(21, 7);
    for r = 1:21
        electrode_grid_map(r, 7) = 100 + r;          % strip 1 (posterior), col A
        electrode_grid_map(r, 6) = 100 + 21 + r;     % strip 1 (posterior), col B
        electrode_grid_map(r, 5) = 100 + 42 + r;     % strip 1 (posterior), col C
        electrode_grid_map(r, 3) = 200 + r;           % strip 2 (anterior),  col A
        electrode_grid_map(r, 2) = 200 + 21 + r;      % strip 2 (anterior),  col B
        electrode_grid_map(r, 1) = 200 + 42 + r;      % strip 2 (anterior),  col C
    end

    % Build the UI cell for every electrode (nav arrow + mark button)
    elec_ui = struct();   % keyed by 'e<label>', e.g. elec_ui.e101
    for r = 1:21
        for c = [1 2 3 5 6 7]
            elec_lbl   = electrode_grid_map(r, c);
            tbl_idx    = find(electrodes.name == elec_lbl, 1);
            full_color = electrodes.color(tbl_idx, :);            % assigned jet colour
            dim_color  = full_color * 0.2 + [0.75 0.75 0.75];   % muted colour for unset state

            % Outer container panel (handles panel-background clicks)
            cell_panel = uipanel(strip_grid, 'BorderType','none');
            cell_panel.Layout.Row    = r + 1;   % +1 offset because row 1 is the header
            cell_panel.Layout.Column = c;
            cell_panel.ButtonDownFcn = @(~,~) handle_electrode_panel_click(elec_lbl);

            % Inner grid: [20 px nav arrow | flex mark button]
            inner_grid = uigridlayout(cell_panel, [1 2]);
            inner_grid.ColumnWidth     = {20, '1x'};
            inner_grid.Padding         = [1 1 1 1];
            inner_grid.BackgroundColor = dim_color;

            elec_ui.(sprintf('e%d', elec_lbl)).container  = inner_grid;
            elec_ui.(sprintf('e%d', elec_lbl)).true_color = full_color;

            % Navigation arrow: jumps all three views to this electrode's location
            elec_ui.(sprintf('e%d', elec_lbl)).nav_arrow = uibutton(inner_grid, ...
                'Text','▶', 'BackgroundColor','k', 'FontColor',[0.5 0.5 0.5], 'FontSize',10, ...
                'Tooltip','Go to Electrode', 'Enable','off', ...
                'ButtonPushedFcn', @(~,~) go_to_electrode(elec_lbl));

            % Mark button: shows label number; clicking enters mark-placement mode
            elec_ui.(sprintf('e%d', elec_lbl)).mark_btn = uibutton(inner_grid, ...
                'Text', num2str(elec_lbl), ...
                'BackgroundColor','k', 'FontColor','w', 'FontSize',9, ...
                'Tooltip','Mark Electrode', ...
                'ButtonPushedFcn', @(~,~) activate_mark_mode(elec_lbl));
        end
    end

    % ---- Bottom row: visibility toggles & delete button ----
    bottom_ctrl_row = uigridlayout(electrode_panel_layout, [1 3]);
    bottom_ctrl_row.ColumnWidth = {'1x', '1x', 'fit'};
    bottom_ctrl_row.Padding     = [0 0 0 0];

    cb_out_of_plane = uicheckbox(bottom_ctrl_row, ...
        'Text','Show out-of-plane electrodes (o)', 'Value',1, ...
        'ValueChangedFcn', @(src,~) toggle_out_of_plane(src.Value));

    cb_crosshairs = uicheckbox(bottom_ctrl_row, ...
        'Text','Show view crosshairs (c)', 'Value',1, ...
        'ValueChangedFcn', @(src,~) toggle_crosshairs(src.Value));

    btn_delete_mode = uibutton(bottom_ctrl_row, 'Text','Delete marker', ...
        'ButtonPushedFcn', @(~,~) toggle_delete_mode());

    %% ======================================================================
    %% 6.  RIGHT PANEL — 2x2 BRAIN VIEWER GRID
    %% ======================================================================

    brain_viewer_panel = uipanel(top_layout, 'Title','Subject: [Not Loaded]');
    view_grid = uigridlayout(brain_viewer_panel, [2 2]);

    % Grid cell occupancy (filled left-to-right then top-to-bottom):
    %   [1,1] Coronal (Y)   [1,2] Axial (Z)
    %   [2,1] Sagittal (X)  [2,2] Controls
    [ax_coronal,  sl_coronal,  btn_max_cor] = create_view_panel(view_grid, 'Coronal (Y)');
    [ax_axial,    sl_axial,    btn_max_axi] = create_view_panel(view_grid, 'Axial (Z)');
    [ax_sagittal, sl_sagittal, btn_max_sag] = create_view_panel(view_grid, 'Sagittal (X)');

    % Controls / legend panel (bottom-right quadrant).
    % The textarea width is set generously so the longest line fits without wrapping.
    controls_panel = uipanel(view_grid, 'Title','Controls');
    controls_text  = uitextarea(controls_panel, ...
        'Position',[10 10 420 200], 'Editable','off', ...
        'BackgroundColor',[0 0 0], 'FontColor',[1 1 1]);
    controls_text.Value = { ...
        'Mouse wheel scroll: shift view plane', ...
        'Shift + mouse wheel scroll: shift view plane (4x speed)', ...
        'e: Mark next electrode', ...
        'd: Delete marker mode', ...
        'o: Toggle out-of-plane electrode visibility', ...
        'c: Toggle view crosshairs', ...
        'm: Maximize/restore hovered view', ...
        'Ctrl+o: Load electrode locations', ...
        'Ctrl+s: Save electrode locations'};

    %% ======================================================================
    %% NESTED FUNCTION DEFINITIONS
    %% ======================================================================

    % ---- File I/O ----------------------------------------------------------

    function load_brain_volume()
        % Try to load from the expected BIDS-style path; fall back to a file picker.
        default_path = fullfile(options.root_dir, 'DBS', 'derivatives', ...
            ['sub-' options.subject_id], 'ecogloc', 'rpreop', 'rpreop_ct.nii');

        if ~exist(default_path, 'file')
            [chosen_file, chosen_dir] = uigetfile('*.nii', 'Select CT Volume', options.root_dir);
            if chosen_file == 0, return; end    % user cancelled
            default_path = fullfile(chosen_dir, chosen_file);
        end

        try
            prog = uiprogressdlg(main_fig, 'Title','Loading', 'Indeterminate','on', ...
                'Message', sprintf('Loading CT volume:\n%s', default_path));
            drawnow;

            nifti_info   = niftiinfo(default_path);
            brain_volume = double(niftiread(nifti_info));  % load as double for imagesc scaling
            volume_size  = size(brain_volume);             % [nX, nY, nZ]
            current_slices = round(volume_size / 2);       % start viewing at the centre

            % Configure each slice slider for the loaded volume dimensions
            configure_slice_slider(sl_sagittal, volume_size(1), current_slices(1));
            configure_slice_slider(sl_coronal,  volume_size(2), current_slices(2));
            configure_slice_slider(sl_axial,    volume_size(3), current_slices(3));

            brain_viewer_panel.Title = sprintf('Subject: %s     File: %s', ...
                options.subject_id, default_path);
            update_views();
            close(prog);
        catch err
            if exist('prog','var'), close(prog); end
            uialert(main_fig, ['Failed to load NIfTI: ' err.message], 'Load Error');
        end
    end

    function load_electrode_locations()
        % Load previously-saved electrode coordinates from a .tsv file.
        if isempty(brain_volume)
            uialert(main_fig, ...
                'Please load the brain volume before loading electrode locations.', 'Warning');
            return;
        end
        tsv_dir = fullfile(options.root_dir, 'DBS', 'derivatives', ...
            ['sub-' options.subject_id], 'ecogloc');
        [chosen_file, chosen_dir] = uigetfile('*.tsv', 'Load TSV', tsv_dir);
        if chosen_file == 0, return; end    % user cancelled

        try
            loaded = readtable(fullfile(chosen_dir, chosen_file), ...
                'FileType','text', 'Delimiter','\t');

            electrodes.voxel_x = loaded.voxel_x;
            electrodes.voxel_y = loaded.voxel_y;
            electrodes.voxel_z = loaded.voxel_z;

            if ismember('mm_x', loaded.Properties.VariableNames)
                % File already contains mm coords — use them directly
                electrodes.mm_x = loaded.mm_x;
                electrodes.mm_y = loaded.mm_y;
                electrodes.mm_z = loaded.mm_z;
            else
                % Recompute mm coords from voxel coords using the NIfTI affine transform
                for row_i = 1:height(loaded)
                    if ~isnan(loaded.voxel_x(row_i)) && ...
                       isfield(nifti_info,'Transform') && ~isempty(nifti_info.Transform)
                        [electrodes.mm_x(row_i), electrodes.mm_y(row_i), electrodes.mm_z(row_i)] = ...
                            transformPointsForward(nifti_info.Transform, ...
                                loaded.voxel_x(row_i), loaded.voxel_y(row_i), loaded.voxel_z(row_i));
                    end
                end
            end

            % Restore button text (checkmark) and container colour for placed electrodes
            for ei = 1:num_electrodes
                lbl = electrodes.name(ei);
                ui  = elec_ui.(sprintf('e%d', lbl));
                if ~isnan(electrodes.voxel_x(ei))
                    ui.mark_btn.Text             = [num2str(lbl) char(10003)];  % add checkmark
                    ui.container.BackgroundColor = ui.true_color;
                else
                    ui.mark_btn.Text             = num2str(lbl);
                    ui.container.BackgroundColor = ui.true_color * 0.2 + [0.75 0.75 0.75];
                end
            end

            unsaved_electrode_changes = false;
            update_views();
        catch err
            uialert(main_fig, ['Failed to load TSV: ' err.message], 'Error');
        end
        focus_figure();
    end

    function save_electrode_locations()
        % Save all electrode coordinates to a .tsv file.
        % Numeric labels are exported as prefixed strings (e.g. 101 -> 'ecog_L101').
        default_filename = sprintf('sub-%s_eleclocalizer-voxel-coords.tsv', options.subject_id);
        save_path = fullfile(options.root_dir, 'DBS', 'derivatives', ...
            ['sub-' options.subject_id], 'ecogloc', default_filename);
        [chosen_file, chosen_dir] = uiputfile('*.tsv', 'Save TSV', save_path);
        if chosen_file == 0, return; end    % user cancelled

        % Work on a copy so the in-memory table (which uses numeric labels) is unchanged
        export_tbl      = electrodes;
        export_tbl.name = arrayfun(@(x) sprintf('ecog_L%d', x), export_tbl.name, ...
            'UniformOutput', false);

        writetable(export_tbl, fullfile(chosen_dir, chosen_file), ...
            'FileType','text', 'Delimiter','\t');
        unsaved_electrode_changes = false;
        focus_figure();
    end

    % ---- View Creation -----------------------------------------------------

    function [axes_out, slider_out, btn_maximize] = create_view_panel(parent_grid, view_title)
        % create_view_panel  Build one orthogonal brain-view panel inside parent_grid.
        %
        % Panel is a 2-row x 2-column grid:
        %   [Row 1, Col 1]    top bar   — title label | Flip up/down toggle | Flip R/L toggle | Maximize btn
        %   [Row 2, Col 1]    uiaxes    — displays the 2-D brain slice image
        %   [Rows 1-2, Col 2] slider column — [▲ | vertical slice slider | ▼]
        %
        % The two flip buttons are state (toggle) buttons.  Their ValueChangedFcn callbacks
        % write into the shared flip_*_view_ud / flip_*_view_lr flags and call update_views().
        % They are NOT returned as outputs because all their state is managed through the
        % shared nested-function workspace; no external code needs to reference them.
        %
        % Returns:
        %   axes_out     – uiaxes handle for the brain slice image
        %   slider_out   – vertical uislider for navigating through slices
        %   btn_maximize – push button to maximise / restore this view
        %                  (returned so toggle_maximize can update its label text)

        panel_grid = uigridlayout(parent_grid, [2 2]);
        panel_grid.RowHeight   = {25, '1x'};
        panel_grid.ColumnWidth = {'1x', 25};

        % Top bar: [title label (flex) | Flip up/down (100px) | Flip R/L (80px) | Maximize (110px)]
        top_bar = uigridlayout(panel_grid, [1 4]);
        top_bar.Layout.Row    = 1;
        top_bar.Layout.Column = 1;
        top_bar.ColumnWidth   = {'1x', 100, 80, 110};
        top_bar.Padding       = [0 0 0 0];

        uilabel(top_bar, 'Text', view_title, 'FontWeight','bold');

        % "Flip up/down": toggles YDir between 'reverse' (default) and 'normal'.
        % Useful when the brain appears vertically inverted due to acquisition orientation.
        uibutton(top_bar, 'state', ...
            'Text','Flip up/down', 'FontSize',10, ...
            'Tooltip',['Flip the vertical display orientation of this view. ' ...
                       'Does not affect stored voxel or mm coordinates.'], ...
            'ValueChangedFcn', @(src,~) toggle_view_flip_ud(view_title, src.Value));

        % "Flip R/L": toggles XDir between 'normal' (default) and 'reverse'.
        % Useful when the brain appears left-right mirrored relative to expectation.
        uibutton(top_bar, 'state', ...
            'Text','Flip R/L', 'FontSize',10, ...
            'Tooltip',['Flip the horizontal display orientation of this view. ' ...
                       'Does not affect stored voxel or mm coordinates.'], ...
            'ValueChangedFcn', @(src,~) toggle_view_flip_lr(view_title, src.Value));

        btn_maximize = uibutton(top_bar, 'Text','Maximize (m)', 'FontSize',11, ...
            'ButtonPushedFcn', @(~,~) toggle_maximize(view_title));

        % Brain image axes
        axes_out = uiaxes(panel_grid);
        axes_out.Layout.Row    = 2;
        axes_out.Layout.Column = 1;
        disableDefaultInteractivity(axes_out);  % suppress zoom / pan / data-tips

        % Vertical slider column (spans both rows, right side of the panel)
        slider_col = uigridlayout(panel_grid, [3 1]);
        slider_col.Layout.Row    = [1 2];   % span both rows
        slider_col.Layout.Column = 2;
        slider_col.RowHeight     = {25, '1x', 25};
        slider_col.Padding       = [0 0 0 0];
        slider_col.RowSpacing    = 2;

        uibutton(slider_col, 'Text','▲', ...
            'ButtonPushedFcn', @(~,~) shift_slice_by_button(view_title, +1));
        slider_out = uislider(slider_col, 'Orientation','vertical', ...
            'ValueChangedFcn', @(src,~) slider_moved(view_title, src.Value));
        uibutton(slider_col, 'Text','▼', ...
            'ButtonPushedFcn', @(~,~) shift_slice_by_button(view_title, -1));
    end

    % ---- Slice Navigation --------------------------------------------------

    function configure_slice_slider(sl, max_slice, init_slice)
        % Set a slice slider's range and initial value; hide tick marks.
        sl.Limits     = [1, max_slice];
        sl.Value      = init_slice;
        sl.MajorTicks = [];
        sl.MinorTicks = [];
    end

    function shift_slice_by_button(view_name, direction)
        % Advance or retreat the slice index by one step via the ▲ / ▼ buttons.
        switch view_name
            case 'Sagittal (X)',  current_slices(1) = current_slices(1) + direction;
            case 'Coronal (Y)',   current_slices(2) = current_slices(2) + direction;
            case 'Axial (Z)',     current_slices(3) = current_slices(3) + direction;
        end
        update_views();
        focus_figure();
    end

    function slider_moved(view_name, new_value)
        % Called when the user drags a slice slider; snap to the nearest integer slice.
        snapped_value = round(new_value);
        switch view_name
            case 'Sagittal (X)',  current_slices(1) = snapped_value;
            case 'Coronal (Y)',   current_slices(2) = snapped_value;
            case 'Axial (Z)',     current_slices(3) = snapped_value;
        end
        update_views();
        focus_figure();
    end

    % ---- Core Display ------------------------------------------------------

    function update_views()
        % Redraw all three brain slice views from scratch.
        % Called whenever the slice position, electrode data, or display options change.
        if isempty(brain_volume), return; end

        % Clamp slice indices to valid volume bounds before rendering
        current_slices = max(1, min(current_slices, volume_size));

        % Keep sliders in sync with current slice indices
        sl_sagittal.Value = current_slices(1);
        sl_coronal.Value  = current_slices(2);
        sl_axial.Value    = current_slices(3);

        % Extract 2-D slice matrices from the 3-D volume.
        % brain_volume dimensions: (X, Y, Z)
        %   Sagittal slice (fix X): rows = Y,  cols = Z  ->  size [nY x nZ]
        %   Coronal  slice (fix Y): rows = X,  cols = Z  ->  size [nX x nZ]
        %   Axial    slice (fix Z): rows = X,  cols = Y  ->  size [nX x nY]
        slice_sag = squeeze(brain_volume(current_slices(1), :, :));  % Y x Z
        slice_cor = squeeze(brain_volume(:, current_slices(2), :));  % X x Z
        slice_axi = squeeze(brain_volume(:, :, current_slices(3)));  % X x Y

        % ---- Sagittal view (perpendicular to X axis) ----
        cla(ax_sagittal);
        im_sag = imagesc(ax_sagittal, slice_sag);
        colormap(ax_sagittal, 'gray');
        % Apply display orientation flags.  These must be set AFTER imagesc because imagesc
        % internally resets YDir to 'reverse' and XDir to 'normal' on every call.
        ax_sagittal.YDir = ifelse(flip_sagittal_view_ud, 'normal',  'reverse');
        ax_sagittal.XDir = ifelse(flip_sagittal_view_lr, 'reverse', 'normal');
        % Disable hit-testing on the image object so clicks pass through to the axes
        im_sag.HitTest       = 'off';
        im_sag.PickableParts = 'none';
        ax_sagittal.ButtonDownFcn = @(~,~) axis_clicked('Sagittal', ax_sagittal.CurrentPoint);
        ax_sagittal.HitTest       = 'on';
        ax_sagittal.PickableParts = 'all';

        % ---- Coronal view (perpendicular to Y axis) ----
        cla(ax_coronal);
        im_cor = imagesc(ax_coronal, slice_cor);
        colormap(ax_coronal, 'gray');
        ax_coronal.YDir = ifelse(flip_coronal_view_ud, 'normal',  'reverse');
        ax_coronal.XDir = ifelse(flip_coronal_view_lr, 'reverse', 'normal');
        im_cor.HitTest       = 'off';
        im_cor.PickableParts = 'none';
        ax_coronal.ButtonDownFcn = @(~,~) axis_clicked('Coronal', ax_coronal.CurrentPoint);
        ax_coronal.HitTest       = 'on';
        ax_coronal.PickableParts = 'all';

        % ---- Axial view (perpendicular to Z axis) ----
        cla(ax_axial);
        im_axi = imagesc(ax_axial, slice_axi);
        colormap(ax_axial, 'gray');
        ax_axial.YDir = ifelse(flip_axial_view_ud, 'normal',  'reverse');
        ax_axial.XDir = ifelse(flip_axial_view_lr, 'reverse', 'normal');
        im_axi.HitTest       = 'off';
        im_axi.PickableParts = 'none';
        ax_axial.ButtonDownFcn = @(~,~) axis_clicked('Axial', ax_axial.CurrentPoint);
        ax_axial.HitTest       = 'on';
        ax_axial.PickableParts = 'all';

        % ---- Crosshairs ----
        % Each view shows two dashed lines at the intersection of the other two slice planes.
        % yline / xline use data coordinates and are therefore unaffected by YDir / XDir.
        if show_crosshairs
            hold(ax_sagittal,'on');
            % Sagittal axes: rows = Y, cols = Z
            yline(ax_sagittal, current_slices(2), '--', 'Color',[0.7 0.7 0.7]);  % coronal  (Y) plane
            xline(ax_sagittal, current_slices(3), '--', 'Color',[0.7 0.7 0.7]);  % axial    (Z) plane
            hold(ax_sagittal,'off');

            hold(ax_coronal,'on');
            % Coronal axes: rows = X, cols = Z
            yline(ax_coronal, current_slices(1), '--', 'Color',[0.7 0.7 0.7]);   % sagittal (X) plane
            xline(ax_coronal, current_slices(3), '--', 'Color',[0.7 0.7 0.7]);   % axial    (Z) plane
            hold(ax_coronal,'off');

            hold(ax_axial,'on');
            % Axial axes: rows = X, cols = Y
            yline(ax_axial, current_slices(1), '--', 'Color',[0.7 0.7 0.7]);     % sagittal (X) plane
            xline(ax_axial, current_slices(2), '--', 'Color',[0.7 0.7 0.7]);     % coronal  (Y) plane
            hold(ax_axial,'off');
        end

        draw_electrode_markers();
        update_arrow_states();
    end

    function draw_electrode_markers()
        % Overlay scatter-plot markers for all placed electrodes on the three views.
        %
        % Marker-size convention:
        %   60 pt^2 = electrode is IN the current slice plane  (distance < 1 voxel)
        %   20 pt^2 = electrode is out of plane (only drawn when show_out_of_plane = true)
        %
        % Coordinate mapping per view (imagesc rows/cols serve as axis data coords).
        % These mappings are unaffected by XDir / YDir; MATLAB's CurrentPoint and scatter()
        % both work in the same data-coordinate space regardless of axis direction settings.
        %   Sagittal  (rows=Y, cols=Z) -> scatter(ax, vox_z, vox_y)
        %   Coronal   (rows=X, cols=Z) -> scatter(ax, vox_z, vox_x)
        %   Axial     (rows=X, cols=Y) -> scatter(ax, vox_y, vox_x)

        hold(ax_sagittal,'on');
        hold(ax_coronal, 'on');
        hold(ax_axial,   'on');

        for ei = 1:num_electrodes
            if isnan(electrodes.voxel_x(ei)), continue; end  % electrode not yet placed

            vox_x = electrodes.voxel_x(ei);
            vox_y = electrodes.voxel_y(ei);
            vox_z = electrodes.voxel_z(ei);
            ec    = electrodes.color(ei, :);

            % Distance (voxels) from this electrode to each current slice plane
            dist_from_sag_plane = abs(vox_x - current_slices(1));
            dist_from_cor_plane = abs(vox_y - current_slices(2));
            dist_from_axi_plane = abs(vox_z - current_slices(3));

            % Sagittal view
            if dist_from_sag_plane < 1 || show_out_of_plane
                msz = ifelse(dist_from_sag_plane < 1, 60, 20);
                s = scatter(ax_sagittal, vox_z, vox_y, msz, ec, 'filled','MarkerEdgeColor','k');
                s.HitTest       = 'on';
                s.PickableParts = 'all';
                s.ButtonDownFcn = @(~,~) axis_clicked('Sagittal', ax_sagittal.CurrentPoint);
            end

            % Coronal view
            if dist_from_cor_plane < 1 || show_out_of_plane
                msz = ifelse(dist_from_cor_plane < 1, 60, 20);
                s = scatter(ax_coronal, vox_z, vox_x, msz, ec, 'filled','MarkerEdgeColor','k');
                s.HitTest       = 'on';
                s.PickableParts = 'all';
                s.ButtonDownFcn = @(~,~) axis_clicked('Coronal', ax_coronal.CurrentPoint);
            end

            % Axial view
            if dist_from_axi_plane < 1 || show_out_of_plane
                msz = ifelse(dist_from_axi_plane < 1, 60, 20);
                s = scatter(ax_axial, vox_y, vox_x, msz, ec, 'filled','MarkerEdgeColor','k');
                s.HitTest       = 'on';
                s.PickableParts = 'all';
                s.ButtonDownFcn = @(~,~) axis_clicked('Axial', ax_axial.CurrentPoint);
            end
        end

        hold(ax_sagittal,'off');
        hold(ax_coronal, 'off');
        hold(ax_axial,   'off');
    end

    function update_arrow_states()
        % Refresh the navigation arrow (▶) appearance for every electrode.
        %
        % Visual states:
        %   Disabled / grey      -> electrode has not been placed
        %   Large white-on-white -> electrode is in (or very near) the current slice plane
        %   Small white-on-black -> electrode is placed but not in any current slice plane

        for ei = 1:num_electrodes
            lbl = electrodes.name(ei);
            ui  = elec_ui.(sprintf('e%d', lbl));

            if isnan(electrodes.voxel_x(ei))
                % Not yet placed
                ui.nav_arrow.FontColor       = [0.5 0.5 0.5];
                ui.nav_arrow.BackgroundColor = 'k';
                ui.nav_arrow.FontSize        = 10;
                ui.nav_arrow.Enable          = 'off';
            else
                % Placed: enable arrow and show electrode's full colour in the container
                ui.nav_arrow.Enable          = 'on';
                ui.container.BackgroundColor = ui.true_color;

                % Highlight if visible in at least one current slice plane
                in_current_plane = ...
                    abs(electrodes.voxel_x(ei) - current_slices(1)) < 1 || ...
                    abs(electrodes.voxel_y(ei) - current_slices(2)) < 1 || ...
                    abs(electrodes.voxel_z(ei) - current_slices(3)) < 1;

                if in_current_plane
                    ui.nav_arrow.FontColor       = 'k';
                    ui.nav_arrow.BackgroundColor = 'w';
                    ui.nav_arrow.FontSize        = 16;
                else
                    ui.nav_arrow.FontColor       = 'w';
                    ui.nav_arrow.BackgroundColor = 'k';
                    ui.nav_arrow.FontSize        = 10;
                end
            end
        end
    end

    % ---- Electrode Interaction ---------------------------------------------

    function activate_mark_mode(elec_lbl)
        % Enter placement mode for the specified electrode.
        % The next click in any brain view will assign that electrode's location.
        % If delete mode is currently active, deletes the electrode instead.
        if delete_mode
            delete_electrode(elec_lbl);
            return;
        end

        % Deactivate the previously active electrode (restore its dim background)
        if ~isnan(active_electrode)
            prev_ui = elec_ui.(sprintf('e%d', active_electrode));
            prev_ui.container.BackgroundColor = prev_ui.true_color * 0.2 + [0.75 0.75 0.75];
        end

        active_electrode = elec_lbl;
        main_fig.Pointer = 'crosshair';   % visual cue: waiting for a placement click

        new_ui = elec_ui.(sprintf('e%d', elec_lbl));
        new_ui.container.BackgroundColor = 'w';   % highlight the newly active electrode
        focus_figure();
    end

    function axis_clicked(view_name, cursor_pt)
        % Handle a mouse click inside one of the brain-view axes.
        %
        % DELETE mode -> find the placed electrode closest to the click and remove it.
        % MARK mode   -> place the active electrode at the clicked voxel position.
        %
        % cursor_pt is ax.CurrentPoint (size [2x3]); we use row 1, columns 1-2:
        %   cursor_pt(1,1) = horizontal data coordinate (column index in the slice image)
        %   cursor_pt(1,2) = vertical   data coordinate (row    index in the slice image)
        %
        % Crucially, these data coordinates are always in the image's row/column index space
        % regardless of YDir or XDir, so neither the up/down nor the R/L flip toggles require
        % any adjustment here — clicks always map to the correct voxel automatically.
        %
        % View-to-voxel mapping:
        %   Sagittal (rows=Y, cols=Z, fix X): col->vox_z, row->vox_y, vox_x = current_slices(1)
        %   Coronal  (rows=X, cols=Z, fix Y): col->vox_z, row->vox_x, vox_y = current_slices(2)
        %   Axial    (rows=X, cols=Y, fix Z): col->vox_y, row->vox_x, vox_z = current_slices(3)

        if delete_mode
            % Reconstruct the approximate 3-D voxel position of the click
            switch view_name
                case 'Sagittal'
                    click_x = current_slices(1);
                    click_y = round(cursor_pt(1,2));    % row  -> Y
                    click_z = round(cursor_pt(1,1));    % col  -> Z
                case 'Coronal'
                    click_y = current_slices(2);
                    click_x = round(cursor_pt(1,2));    % row  -> X
                    click_z = round(cursor_pt(1,1));    % col  -> Z
                case 'Axial'
                    click_z = current_slices(3);
                    click_x = round(cursor_pt(1,2));    % row  -> X
                    click_y = round(cursor_pt(1,1));    % col  -> Y
                otherwise, return;
            end

            % Find the nearest placed electrode to the click position
            min_dist      = inf;
            nearest_label = NaN;
            for ei = 1:num_electrodes
                if isnan(electrodes.voxel_x(ei)), continue; end
                d = sqrt((electrodes.voxel_x(ei) - click_x)^2 + ...
                         (electrodes.voxel_y(ei) - click_y)^2 + ...
                         (electrodes.voxel_z(ei) - click_z)^2);
                if d < min_dist
                    min_dist      = d;
                    nearest_label = electrodes.name(ei);
                end
            end
            % Only delete if the click was within 10 voxels of an electrode
            if min_dist <= 10
                delete_electrode(nearest_label);
            end
            return;
        end

        % MARK mode: place the active electrode at the clicked location
        if isnan(active_electrode), return; end

        col_val = round(cursor_pt(1,1));   % horizontal axis value
        row_val = round(cursor_pt(1,2));   % vertical   axis value

        switch view_name
            case 'Sagittal'
                new_vox_x = current_slices(1);  new_vox_y = row_val;  new_vox_z = col_val;
            case 'Coronal'
                new_vox_y = current_slices(2);  new_vox_x = row_val;  new_vox_z = col_val;
            case 'Axial'
                new_vox_z = current_slices(3);  new_vox_x = row_val;  new_vox_y = col_val;
            otherwise, return;
        end

        % Store voxel coordinates in the master table
        ti = find(electrodes.name == active_electrode, 1);
        electrodes.voxel_x(ti) = new_vox_x;
        electrodes.voxel_y(ti) = new_vox_y;
        electrodes.voxel_z(ti) = new_vox_z;

        % Compute millimetre coordinates using the NIfTI affine transform
        if isfield(nifti_info,'Transform') && ~isempty(nifti_info.Transform)
            [mm_x, mm_y, mm_z] = transformPointsForward( ...
                nifti_info.Transform, new_vox_x, new_vox_y, new_vox_z);
        elseif isfield(nifti_info,'PixelDimensions')
            % Fallback when no rotation is present: scale by voxel size
            mm_x = new_vox_x * nifti_info.PixelDimensions(1);
            mm_y = new_vox_y * nifti_info.PixelDimensions(2);
            mm_z = new_vox_z * nifti_info.PixelDimensions(3);
        else
            mm_x = NaN; mm_y = NaN; mm_z = NaN;
        end
        electrodes.mm_x(ti) = mm_x;
        electrodes.mm_y(ti) = mm_y;
        electrodes.mm_z(ti) = mm_z;

        % Update the electrode's mark button to display a checkmark
        ui = elec_ui.(sprintf('e%d', active_electrode));
        ui.mark_btn.Text             = [num2str(active_electrode) char(10003)];
        ui.container.BackgroundColor = ui.true_color;

        unsaved_electrode_changes = true;
        active_electrode = NaN;         % exit mark mode
        main_fig.Pointer = 'arrow';
        update_views();
    end

    function delete_electrode(elec_lbl)
        % Clear all stored coordinate data for one electrode and reset its UI state.
        ti = find(electrodes.name == elec_lbl, 1);
        if isempty(ti) || isnan(electrodes.voxel_x(ti)), return; end  % nothing to delete

        % Clear both voxel and mm coordinates
        electrodes.voxel_x(ti) = NaN;  electrodes.mm_x(ti) = NaN;
        electrodes.voxel_y(ti) = NaN;  electrodes.mm_y(ti) = NaN;
        electrodes.voxel_z(ti) = NaN;  electrodes.mm_z(ti) = NaN;

        ui = elec_ui.(sprintf('e%d', elec_lbl));
        ui.mark_btn.Text             = num2str(elec_lbl);                         % remove checkmark
        ui.container.BackgroundColor = ui.true_color * 0.2 + [0.75 0.75 0.75];   % dim background
        ui.nav_arrow.Enable          = 'off';

        unsaved_electrode_changes = true;
        update_views();
    end

    function go_to_electrode(elec_lbl)
        % Navigate all three slice views to the planes that contain this electrode.
        ti = find(electrodes.name == elec_lbl, 1);
        if ~isempty(ti) && ~isnan(electrodes.voxel_x(ti))
            current_slices = round([electrodes.voxel_x(ti), ...
                                    electrodes.voxel_y(ti), ...
                                    electrodes.voxel_z(ti)]);
            update_views();
        end
        focus_figure();
    end

    function handle_electrode_panel_click(elec_lbl)
        % Handle a click on the background area of an electrode's panel cell.
        if delete_mode
            delete_electrode(elec_lbl);
        else
            activate_mark_mode(elec_lbl);
        end
    end

    % ---- View Orientation & Layout -----------------------------------------

    function toggle_view_flip_ud(view_name, new_flip_state)
        % Toggle the vertical display orientation for one of the three brain views.
        %
        % ON  (new_flip_state = true)  -> YDir = 'normal'  (row 1 at bottom; image flipped up/down)
        % OFF (new_flip_state = false) -> YDir = 'reverse' (row 1 at top;    imagesc default)
        %
        % This is purely cosmetic — MATLAB's CurrentPoint and scatter() both operate in the
        % axes data-coordinate space regardless of YDir, so voxel placement is unaffected.
        switch view_name
            case 'Sagittal (X)', flip_sagittal_view_ud = new_flip_state;
            case 'Coronal (Y)',  flip_coronal_view_ud  = new_flip_state;
            case 'Axial (Z)',    flip_axial_view_ud    = new_flip_state;
        end
        update_views();
        focus_figure();
    end

    function toggle_view_flip_lr(view_name, new_flip_state)
        % Toggle the horizontal display orientation for one of the three brain views.
        %
        % ON  (new_flip_state = true)  -> XDir = 'reverse' (col 1 at right; image mirrored L/R)
        % OFF (new_flip_state = false) -> XDir = 'normal'  (col 1 at left;  imagesc default)
        %
        % This is purely cosmetic — MATLAB's CurrentPoint and scatter() both operate in the
        % axes data-coordinate space regardless of XDir, so voxel placement is unaffected.
        switch view_name
            case 'Sagittal (X)', flip_sagittal_view_lr = new_flip_state;
            case 'Coronal (Y)',  flip_coronal_view_lr  = new_flip_state;
            case 'Axial (Z)',    flip_axial_view_lr    = new_flip_state;
        end
        update_views();
        focus_figure();
    end

    function toggle_maximize(view_name)
        % Expand the chosen view to fill the full 2x2 brain-viewer grid,
        % or restore all four quadrants to equal size if already maximised.
        if isempty(maximized_view)
            maximized_view = view_name;
            % Hide the three non-selected quadrants by collapsing their row/column to zero
            switch view_name
                case 'Coronal (Y)'     % grid cell [1,1]
                    view_grid.RowHeight   = {'1x', 0};
                    view_grid.ColumnWidth = {'1x', 0};
                case 'Axial (Z)'       % grid cell [1,2]
                    view_grid.RowHeight   = {'1x', 0};
                    view_grid.ColumnWidth = {0, '1x'};
                case 'Sagittal (X)'    % grid cell [2,1]
                    view_grid.RowHeight   = {0, '1x'};
                    view_grid.ColumnWidth = {'1x', 0};
            end
            btn_max_cor.Text = 'Restore (m)';
            btn_max_axi.Text = 'Restore (m)';
            btn_max_sag.Text = 'Restore (m)';
        else
            % Restore equal quadrant sizes
            maximized_view = '';
            view_grid.RowHeight   = {'1x', '1x'};
            view_grid.ColumnWidth = {'1x', '1x'};
            btn_max_cor.Text = 'Maximize (m)';
            btn_max_axi.Text = 'Maximize (m)';
            btn_max_sag.Text = 'Maximize (m)';
        end
        focus_figure();
    end

    % ---- Toggle Helpers ----------------------------------------------------

    function toggle_out_of_plane(new_val)
        show_out_of_plane = new_val;
        update_views();
        focus_figure();
    end

    function toggle_crosshairs(new_val)
        show_crosshairs = new_val;
        update_views();
        focus_figure();
    end

    function toggle_delete_mode()
        % Switch delete mode on or off.
        % While active, the controls panel background turns red as a visual warning,
        % and clicking near any placed electrode removes it.
        delete_mode = ~delete_mode;
        if delete_mode
            main_fig.Pointer                = 'cross';
            controls_text.BackgroundColor   = [0.3 0 0];      % red tint = danger
            btn_delete_mode.BackgroundColor = [0.8 0 0];
            btn_delete_mode.FontColor       = 'w';
            % Cancel any pending electrode placement to avoid conflicting state
            if ~isnan(active_electrode)
                prev_ui = elec_ui.(sprintf('e%d', active_electrode));
                prev_ui.container.BackgroundColor = prev_ui.true_color * 0.2 + [0.75 0.75 0.75];
                active_electrode = NaN;
            end
        else
            main_fig.Pointer                = 'arrow';
            controls_text.BackgroundColor   = [0 0 0];
            btn_delete_mode.BackgroundColor = [0.96 0.96 0.96];
            btn_delete_mode.FontColor       = 'k';
        end
        focus_figure();
    end

    % ---- Event Handlers ----------------------------------------------------

    function handle_keypress(~, key_evt)
        % Keyboard shortcut dispatcher.
        has_ctrl = ismember('control', key_evt.Modifier);
        switch lower(key_evt.Key)
            case 'o'
                if has_ctrl
                    load_electrode_locations();
                else
                    cb_out_of_plane.Value = ~cb_out_of_plane.Value;
                    toggle_out_of_plane(cb_out_of_plane.Value);
                end
            case 's'
                if has_ctrl, save_electrode_locations(); end
            case 'e'
                if ~has_ctrl
                    % Activate mark mode for the first electrode not yet placed
                    next_idx = find(isnan(electrodes.voxel_x), 1, 'first');
                    if ~isempty(next_idx)
                        activate_mark_mode(electrodes.name(next_idx));
                    end
                end
            case 'c'
                if ~has_ctrl
                    cb_crosshairs.Value = ~cb_crosshairs.Value;
                    toggle_crosshairs(cb_crosshairs.Value);
                end
            case 'd'
                if ~has_ctrl, toggle_delete_mode(); end
            case 'm'
                if ~has_ctrl && ~isempty(hovered_view_name)
                    toggle_maximize(hovered_view_name);
                end
        end
    end

    function handle_scroll(~, scroll_evt)
        % Mouse-wheel scroll: advance or retreat the slice in whichever view is under the cursor.
        %
        % Normal scroll (no modifier):
        %   Moves 1 slice per scroll tick  (VerticalScrollCount).
        %
        % Shift + scroll:
        %   Moves 4 slices per scroll tick — useful for coarse navigation through large volumes.
        %   Shift state is read from main_fig.CurrentModifier, which reflects the live
        %   keyboard state at the time of the scroll event.

        scroll_ticks = scroll_evt.VerticalScrollCount;  % positive = scroll down, negative = up

        % Determine scroll speed multiplier based on whether Shift is held
        shift_is_held  = ismember('shift', main_fig.CurrentModifier);
        slices_per_tick = ifelse(shift_is_held, 4, 1);

        % Total slice delta (negative scroll_ticks = wheel up = advance to higher index)
        slice_delta = scroll_ticks * slices_per_tick;

        try
            pt = main_fig.CurrentPoint;  % cursor position in figure pixels

            pos_sag = getpixelposition(ax_sagittal.Parent, true);
            if in_box(pt, pos_sag)
                current_slices(1) = current_slices(1) - slice_delta;
                update_views(); return;
            end

            pos_cor = getpixelposition(ax_coronal.Parent, true);
            if in_box(pt, pos_cor)
                current_slices(2) = current_slices(2) - slice_delta;
                update_views(); return;
            end

            pos_axi = getpixelposition(ax_axial.Parent, true);
            if in_box(pt, pos_axi)
                current_slices(3) = current_slices(3) - slice_delta;
                update_views(); return;
            end
        catch
            % Silently ignore transient errors (e.g. before the first render completes)
        end
    end

    function track_mouse_hover(~, ~)
        % Track which brain view the cursor is currently over.
        %   (1) Updates hovered_view_name, used by the 'm' key shortcut.
        %   (2) In delete mode, draws a hover-ring highlight around the nearest electrode.

        cursor_pos = main_fig.CurrentPoint;
        if isempty(brain_volume), return; end

        hovered_view_name = '';
        hovered_ax        = [];

        pos_sag = getpixelposition(ax_sagittal.Parent, true);
        pos_cor = getpixelposition(ax_coronal.Parent,  true);
        pos_axi = getpixelposition(ax_axial.Parent,    true);

        if in_box(cursor_pos, pos_sag)
            hovered_view_name = 'Sagittal (X)';  hovered_ax = ax_sagittal;
        elseif in_box(cursor_pos, pos_cor)
            hovered_view_name = 'Coronal (Y)';   hovered_ax = ax_coronal;
        elseif in_box(cursor_pos, pos_axi)
            hovered_view_name = 'Axial (Z)';     hovered_ax = ax_axial;
        end

        if delete_mode && ~isempty(hovered_ax)
            % In delete mode: find the nearest electrode to the cursor in the hovered view
            % and draw a ring around it as a visual indicator of what would be deleted.
            ax_cursor = hovered_ax.CurrentPoint(1, 1:2);  % cursor in axes data coordinates

            clear_hover_highlights();
            hovered_electrode = NaN;
            best_dist         = 8;     % hover-detection radius in data (voxel-index) units
            best_screen_pos   = [];

            for ei = 1:num_electrodes
                if isnan(electrodes.voxel_x(ei)), continue; end

                vox_x = electrodes.voxel_x(ei);
                vox_y = electrodes.voxel_y(ei);
                vox_z = electrodes.voxel_z(ei);

                % Convert electrode 3-D voxel coords to the 2-D screen coords of this view,
                % matching the coordinate mapping used in draw_electrode_markers.
                % These are data-space coordinates; XDir/YDir do not affect them.
                if hovered_ax == ax_sagittal
                    dist_from_plane = abs(vox_x - current_slices(1));
                    if dist_from_plane > 5 && ~show_out_of_plane, continue; end
                    screen_x = vox_z;   screen_y = vox_y;   % sagittal: col=Z, row=Y
                elseif hovered_ax == ax_coronal
                    dist_from_plane = abs(vox_y - current_slices(2));
                    if dist_from_plane > 5 && ~show_out_of_plane, continue; end
                    screen_x = vox_z;   screen_y = vox_x;   % coronal:  col=Z, row=X
                else   % ax_axial
                    dist_from_plane = abs(vox_z - current_slices(3));
                    if dist_from_plane > 5 && ~show_out_of_plane, continue; end
                    screen_x = vox_y;   screen_y = vox_x;   % axial:    col=Y, row=X
                end

                hover_dist = sqrt((screen_x - ax_cursor(1))^2 + (screen_y - ax_cursor(2))^2);
                if hover_dist < best_dist
                    best_dist       = hover_dist;
                    best_screen_pos = [screen_x, screen_y, electrodes.name(ei)];
                end
            end

            if ~isempty(best_screen_pos)
                hovered_electrode = best_screen_pos(3);
                hold(hovered_ax,'on');
                % Double ring: white outer for contrast, black inner for definition
                h_outer = plot(hovered_ax, best_screen_pos(1), best_screen_pos(2), 'o', ...
                    'MarkerSize',18, 'LineWidth',3, 'Color','w', ...
                    'LineStyle',':', 'MarkerFaceColor','none', ...
                    'HitTest','off', 'PickableParts','none');
                h_inner = plot(hovered_ax, best_screen_pos(1), best_screen_pos(2), 'o', ...
                    'MarkerSize',14, 'LineWidth',2, 'Color','k', ...
                    'LineStyle',':', 'MarkerFaceColor','none', ...
                    'HitTest','off', 'PickableParts','none');
                hold(hovered_ax,'off');

                % Store handles so they can be deleted on the next hover update
                if     hovered_ax == ax_sagittal,  hover_highlight_handles.sag = [h_outer h_inner];
                elseif hovered_ax == ax_coronal,   hover_highlight_handles.cor = [h_outer h_inner];
                else,                               hover_highlight_handles.axi = [h_outer h_inner];
                end
            end

        else
            % Not in delete mode, or cursor outside all views: clear any lingering rings
            clear_hover_highlights();
            hovered_electrode = NaN;
        end
    end

    function handle_fig_click(~, ~)
        % Figure-level click handler.
        % Intentionally empty: all meaningful click interactions are handled by the
        % per-axes ButtonDownFcn callbacks, which have finer-grained position information.
    end

    % ---- Utilities ---------------------------------------------------------

    function clear_hover_highlights()
        % Delete hover-ring graphics from all three views (used in delete mode).
        if ~isempty(hover_highlight_handles.sag)
            delete(hover_highlight_handles.sag);  hover_highlight_handles.sag = [];
        end
        if ~isempty(hover_highlight_handles.cor)
            delete(hover_highlight_handles.cor);  hover_highlight_handles.cor = [];
        end
        if ~isempty(hover_highlight_handles.axi)
            delete(hover_highlight_handles.axi);  hover_highlight_handles.axi = [];
        end
    end

    function tf = in_box(pt, box)
        % Return true if 2-D point pt = [x, y] lies inside bounding box box = [x0, y0, w, h].
        tf = pt(1) > box(1) && pt(1) < box(1) + box(3) && ...
             pt(2) > box(2) && pt(2) < box(2) + box(4);
    end

    function result = ifelse(condition, true_val, false_val)
        % Inline ternary helper: returns true_val when condition is true, else false_val.
        if condition,  result = true_val;
        else,          result = false_val;
        end
    end

    function focus_figure()
        % Return keyboard focus to the main figure after button / slider interactions.
        drawnow;
        try, figure(main_fig); catch, end
    end

end   % mark_ecog_on_ct