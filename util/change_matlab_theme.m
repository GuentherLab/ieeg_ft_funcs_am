%% on AMD Strix laptop, generally fails to know system theme, but still correctly updates when selected

function change_matlab_theme()
% Detect the current underlying Operating System color scheme
try
    % Query the native Settings API engine directly
    sysTheme = settings.matlab.appearance.MATLABTheme.ActiveValue;
catch
    sysTheme = 'Unknown';
end

% Initialize a modal, clean parent dialog window
fig = uifigure('Name', 'Theme Selector', ...
    'Position', [500, 500, 420, 140], ...
    'WindowStyle', 'modal', ...
    'Resize', 'off');
movegui(fig, 'center'); % Snap directly to screen center

% Add a prompt message label
uilabel(fig, ...
    'Text', 'Select a desktop theme for MATLAB:', ...
    'Position', [20, 95, 380, 25], ...
    'FontWeight', 'bold', ...
    'FontSize', 13);

% Initialize structural choice tracking variable
selectedTheme = '';

% 1. Dark Button (Slightly Dark / Bluish Theme)
btnDark = uibutton(fig, 'push', ...
    'Text', 'Dark', ...
    'Position', [20, 35, 110, 40], ...
    'BackgroundColor', [0.15, 0.22, 0.33], ... % Deep slate blue
    'FontColor', [0.95, 0.95, 0.95], ...       % High contrast crisp white
    'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(src, event) assignChoice('Dark'));

% 2. Light Button (Lighter / Warm Soft Amber Theme)
btnLight = uibutton(fig, 'push', ...
    'Text', 'Light', ...
    'Position', [145, 35, 110, 40], ...
    'BackgroundColor', [0.98, 0.94, 0.88], ... % Soft warm ivory/cream
    'FontColor', [0.20, 0.15, 0.10], ...       % Warm deep charcoal text
    'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(src, event) assignChoice('Light'));

% 3. System Button (Dynamic label based on current system state)
systemLabel = sprintf('System (%s)', sysTheme);

% Pick system button colors matching the actual live state
if strcmpi(sysTheme, 'Dark')
    sysBG = [0.20, 0.20, 0.25];
    sysFG = [0.90, 0.90, 0.90];
else
    sysBG = [0.90, 0.90, 0.92];
    sysFG = [0.10, 0.10, 0.15];
end

btnSystem = uibutton(fig, 'push', ...
    'Text', systemLabel, ...
    'Position', [270, 35, 130, 40], ...
    'BackgroundColor', sysBG, ...
    'FontColor', sysFG, ...
    'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(src, event) assignChoice('System'));

% Halt execution thread until user interacts with the UI buttons
uiwait(fig);

% Embedded callback function to catch user clicks
    function assignChoice(choiceValue)
        selectedTheme = choiceValue;
        close(fig); % Resumes execution after uiwait
    end

% Safely catch window termination/dismissals
if isempty(selectedTheme)
    fprintf('Theme change canceled.\n');
    return;
end

% Commit selection immediately using the App Appearance Engine API
s = settings;
s.matlab.appearance.MATLABTheme.PersonalValue = selectedTheme;

% Print clean diagnostic verification log
fprintf('MATLAB theme successfully changed to: %s\n', selectedTheme);
end
