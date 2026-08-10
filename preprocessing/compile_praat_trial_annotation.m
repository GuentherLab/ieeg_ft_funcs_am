function trials_out = compile_praat_trial_annotation(trials_in, praat_files_dir, cfg)
% COMPILE_PRAAT_TRIAL_ANNOTATION Extracts tier annotations from Praat TextGrid files.
% 
%   trials_out = compile_praat_trial_annotation(trials_in, praat_files_dir)
%   trials_out = compile_praat_trial_annotation(trials_in, praat_files_dir, cfg)

    % 1. Initialization and Input Checking
    if nargin < 3
        cfg = struct();
    end
    
    trials_out = trials_in;
    num_trials = height(trials_in);
    
    % Check if cfg.praat_tiers is provided; if not, infer from the first trial
    if ~isfield(cfg, 'praat_tiers') || isempty(cfg.praat_tiers)
        first_trial_id = trials_in.trial_id(1);
        first_file = fullfile(praat_files_dir, sprintf('trial-%d.textgrid', first_trial_id));
        if ~exist(first_file, 'file')
            error('Default tier generation failed: First trial TextGrid file not found (%s).', first_file);
        end
        tg_data = parse_textgrid(first_file);
        
        % Create default cfg.praat_tiers based on the first file
        num_tiers = length(tg_data);
        cfg.praat_tiers = table('Size', [num_tiers, 4], ...
            'VariableTypes', {'string', 'string', 'cell', 'cell'}, ...
            'VariableNames', {'name', 'class', 'n_times_min_max', 'timepoint_names'});
        for i = 1:num_tiers
            cfg.praat_tiers.name(i) = string(tg_data(i).name);
            cfg.praat_tiers.class(i) = "cell"; % Default missing to 'cell'
            cfg.praat_tiers.n_times_min_max{i} = [];
            cfg.praat_tiers.timepoint_names{i} = {};
        end
    end
    
    % Standardize and validate cfg.praat_tiers
    for i = 1:height(cfg.praat_tiers)
        tier_name = char(cfg.praat_tiers.name(i));
        
        % Handle empty class -> default to 'cell'
        if iscell(cfg.praat_tiers.class)
            t_class = cfg.praat_tiers.class{i};
        else
            t_class = char(cfg.praat_tiers.class(i));
        end
        if isempty(t_class) || ismissing(string(t_class))
            t_class = 'cell';
            if iscell(cfg.praat_tiers.class)
                cfg.praat_tiers.class{i} = t_class;
            else
                cfg.praat_tiers.class(i) = string(t_class);
            end
        end
        
        % Validate timing logic and initialize output variables
        if strcmpi(t_class, 'timing')
            tp_names = cfg.praat_tiers.timepoint_names{i};
            n_min_max = cfg.praat_tiers.n_times_min_max{i};
            
            if ~isempty(tp_names)
                if isempty(n_min_max) || length(n_min_max) ~= 2 || n_min_max(1) ~= n_min_max(2) || n_min_max(1) ~= length(tp_names)
                    error('cfg.praat_tiers validation failed for tier "%s": When timepoint_names is provided, n_times_min_max must contain identical min/max values equal to the number of names.', tier_name);
                end
                % Initialize variables for specified timepoints
                for j = 1:length(tp_names)
                    trials_out.(char(tp_names{j})) = nan(num_trials, 1);
                end
            else
                % Initialize single cell array to hold 1xN doubles
                trials_out.(tier_name) = cell(num_trials, 1);
            end
        else
            % Preallocate for normal classes
            switch lower(t_class)
                case 'double'
                    trials_out.(tier_name) = nan(num_trials, 1);
                case 'string'
                    trials_out.(tier_name) = strings(num_trials, 1);
                case 'logical'
                    % Pre-allocate with NaN for empty checks later (requires double temporarily)
                    % or just use a numeric array initialized to NaN, converting at the very end
                    trials_out.(tier_name) = nan(num_trials, 1);
                case 'cell'
                    trials_out.(tier_name) = cell(num_trials, 1);
                otherwise
                    error('Unsupported class "%s" for tier "%s".', t_class, tier_name);
            end
        end
    end

    % 2. Process Each Trial
    dirfiles = {dir(praat_files_dir).name}'; 
    for itrial = 1:num_trials
        t_id = trials_in.trial_id(itrial);
        gtc_start = double(trials_in.praat_file_start(itrial));

        trialfile = dirfiles(endsWith(dirfiles,sprintf('trial-%d.textgrid', t_id),'IgnoreCase', true)); 
        praat_filepath = fullfile(praat_files_dir, trialfile{1});
        if ~exist(praat_filepath, 'file')
            error('TextGrid for trial_id %d not found.', t_id);
            continue;
        end
        
        % Parse TextGrid
        tg_data = parse_textgrid(praat_filepath);
        
        % Match and process specified tiers
        for itier = 1:height(cfg.praat_tiers)
            tier_name = char(cfg.praat_tiers.name(itier));
            if iscell(cfg.praat_tiers.class)
                t_class = char(cfg.praat_tiers.class{itier});
            else
                t_class = char(cfg.praat_tiers.class(itier));
            end
            n_min_max = cfg.praat_tiers.n_times_min_max{itier};
            tp_names = cfg.praat_tiers.timepoint_names{itier};
            
            % Find matching tier in file
            idx = find(strcmp({tg_data.name}, tier_name), 1);
            if isempty(idx)
                continue; 
            end
            
            tier_info = tg_data(idx);
            num_intervals = length(tier_info.intervals);
            
            % Check for empty condition (1 interval and text is empty)
            if num_intervals == 1 && isempty(strtrim(tier_info.intervals(1).text))
                continue; % Leave pre-allocated spot empty (NaN, <missing>, or [])
            end
            
            % Handle data based on specified class
            switch lower(t_class)
                case {'double', 'string', 'logical'}
                    if num_intervals > 1
                        val_str = strjoin({tier_info.intervals.text}, ', ');
                        error('Single-value assumption violated.\nTrial: %d\nTier: %s\nNumber of values found: %d\nValues: %s', ...
                            t_id, tier_name, num_intervals, val_str);
                    end
                    
                    raw_text = tier_info.intervals(1).text;
                    
                    if strcmpi(t_class, 'double')
                        num_val = str2double(raw_text);
                        if ~isnan(num_val)
                            trials_out.(tier_name)(itrial) = num_val;
                        end
                    elseif strcmpi(t_class, 'string')
                        trials_out.(tier_name)(itrial) = string(raw_text);
                    elseif strcmpi(t_class, 'logical')
                        % assuming text indicates true/false or 1/0
                        if any(strcmpi(raw_text, {'true', '1', 'yes', 't'}))
                            trials_out.(tier_name)(itrial) = 1;
                        elseif any(strcmpi(raw_text, {'false', '0', 'no', 'f'}))
                            trials_out.(tier_name)(itrial) = 0;
                        end
                    end
                    
                case 'cell'
                    cell_vals = {tier_info.intervals.text};
                    trials_out.(tier_name){itrial} = cell_vals;
                    
                case 'timing'
                    if num_intervals <= 1
                        error('Timing extraction error. No timepoints manually marked.\nTrial: %d\nTier: %s', t_id, tier_name);
                    end
                    
                    % Number of marked points is intervals minus 1 (using xmax of all but last)
                    num_marked = num_intervals - 1;
                    
                    if ~isempty(n_min_max) && (num_marked < n_min_max(1) || num_marked > n_min_max(2))
                        error('Trial %d: Tier "%s" has %d timepoints, which is outside the specified range of [%d, %d].', ...
                            t_id, tier_name, num_marked, n_min_max(1), n_min_max(2));
                    end
                    
                    % Calculate GTC for each point
                    xmax_vals = [tier_info.intervals(1:end-1).xmax];
                    gtc_vals = xmax_vals + gtc_start;
                    
                    if isempty(tp_names)
                        % Store as 1xN double in cell
                        trials_out.(tier_name){itrial} = gtc_vals;
                    else
                        % Store in individual named columns
                        for j = 1:length(tp_names)
                            trials_out.(char(tp_names{j}))(itrial) = gtc_vals(j);
                        end
                    end
            end
        end
    end
    
    % Post-process logic types to true logical array (keeping NaN as false, or user handled)
    for i = 1:height(cfg.praat_tiers)
        if iscell(cfg.praat_tiers.class)
            t_class = char(cfg.praat_tiers.class{i});
        else
            t_class = char(cfg.praat_tiers.class(i));
        end
        if strcmpi(t_class, 'logical')
            tier_name = char(cfg.praat_tiers.name(i));
            % Convert double array initialized with NaNs to logical. NaN becomes false.
            nan_mask = isnan(trials_out.(tier_name));
            temp = trials_out.(tier_name);
            temp(nan_mask) = 0; % Default missing logic to false
            trials_out.(tier_name) = logical(temp);
        end
    end

end

% --- Helper Function to Parse TextGrid Files ---
function tg = parse_textgrid(filepath)
    % Reads a TextGrid file and parses IntervalTiers into a struct array
    txt = fileread(filepath);
    
    % Regex to find Tier blocks
    tier_expr = 'item \[\d+\]:\s*class = "(.*?)"\s*name = "(.*?)"\s*xmin = [\d\.]+\s*xmax = [\d\.]+\s*intervals: size = (\d+)(.*?)(?=item \[\d+\]:|$)';
    tier_tokens = regexp(txt, tier_expr, 'tokens');
    
    tg = struct('class', {}, 'name', {}, 'intervals', {});
    
    % Regex to find specific intervals inside a Tier block
    int_expr = 'intervals \[\d+\]:\s*xmin = [\d\.]+\s*xmax = ([\d\.]+)\s*text = "(.*?)"';
    
    for i = 1:length(tier_tokens)
        tg(i).class = tier_tokens{i}{1};
        tg(i).name = tier_tokens{i}{2};
        
        int_block = tier_tokens{i}{4};
        int_tokens = regexp(int_block, int_expr, 'tokens');
        
        intervals = struct('xmax', {}, 'text', {});
        for j = 1:length(int_tokens)
            intervals(j).xmax = str2double(int_tokens{j}{1});
            intervals(j).text = int_tokens{j}{2};
        end
        tg(i).intervals = intervals;
    end
end