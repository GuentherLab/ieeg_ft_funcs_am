% process electrodes:
.... get mean response within specified time windows per trial
.... get warped timecourses based on epochs table
.... add anatomical info from electrodes table to resp table (if op.electrodes_file is provided)

function [resp, trials] = get_epoched_responses(D_in, trials, op)
    ntrials = height(trials);
    field_default('op', 'keep_unwarped_timecourse', false)
    field_default('op', 'trials_to_analyze', true(ntrials, 1))
    field_default('op', 'electrodes_file', '')
    
    % Validate op.trials_to_analyze
    assert(length(op.trials_to_analyze) == ntrials, ...
        'op.trials_to_analyze must be a logical array with the same height as trials.');
    op.trials_to_analyze = logical(op.trials_to_analyze(:));

    epochs = op.epochs; 
    nepochs = height(epochs); 
    % Ensure 'base' is present in op.epochs
    assert(ismember('base', epochs.Properties.RowNames), ...
        '''base'' must be a row in op.epochs');
    
    nchans = length(D_in.label);
    nans_tr = nan(ntrials, 1); 
    false_tr = false(ntrials, 1);
    cel_trials = cell(ntrials, 1); 
    cel_chans_trials_nan = repmat({nans_tr}, nchans, 1);
    cel_chans_trials_false = repmat({false_tr}, nchans, 1);
    
    resp = table(D_in.label, cel_chans_trials_nan, repmat({cel_trials}, nchans, 1), repmat({cel_trials}, nchans, 1), cel_chans_trials_false, ....
      'VariableNames', {'chan', 'base', 'timecourse_unwarped', 'timecourse', 'good_trial'}); 


    %% add anatomical info from electrodes table to resp table (if op.electrodes_file is provided)
    if ~isempty(op.electrodes_file)
        electrode_vars_to_copy = {'chan','type','native_x','native_y','native_z',...
            'mni_x','mni_y','mni_z',...
	        'DISTAL_label_1','DISTAL_weight_1','DISTAL_label_2','DISTAL_weight_2','DISTAL_label_3','DISTAL_weight_3',...
            'HCPMMP1_label_1','HCPMMP1_weight_1','HCPMMP1_label_2','HCPMMP1_weight_2','connector'};

        resp = reref_chan_to_electrode_label(resp); % add 'electrode_label' field so we can match to pre-reref elc table
        electrodes = readtable(op.electrodes_file, 'FileType','text', 'Delimiter','tab'); 
        electrode_vars_to_copy = cellstr(electrode_vars_to_copy);
        
        missingVars = setdiff(electrode_vars_to_copy, electrodes.Properties.VariableNames);
        if ~isempty(missingVars)
            warning('electrodes:missingVarToCopy', ...
                'Variables not found in electrodes table and will be skipped: %s', ...
                strjoin(missingVars, ', '));
        end
        validVarsToCopy = setdiff(electrode_vars_to_copy, missingVars, 'stable');
        
        % normalize to string for matching
        elcNames = string(electrodes.name);
        respLabels = string(resp.electrode_label);
        
        % fail fast if any electrode name is duplicated AND used in resp
        [G, uNames] = findgroups(elcNames);
        cnt = splitapply(@numel, elcNames, G);
        dupNames = uNames(cnt > 1);
        if ~isempty(dupNames) && any(ismember(respLabels, dupNames))
            bad = intersect(respLabels, dupNames);
            error('Multiple electrodes rows match label "%s".', bad(1));
        end
        
        % ensure target columns exist in resp
        for j = 1:numel(validVarsToCopy)
            vn = validVarsToCopy{j};
            if ~ismember(vn, resp.Properties.VariableNames)
                template = electrodes.(vn);
                if isnumeric(template) || islogical(template)
                    resp.(vn) = NaN(height(resp), 1);
                elseif isstring(template)
                    resp.(vn) = strings(height(resp),1);
                    resp.(vn)(:) = missing;
                elseif iscell(template)
                    resp.(vn) = repmat({missing}, height(resp),1);
                else
                    resp.(vn) = repmat(missing, height(resp),1);
                end
            end
        end
        
        % copy
        for i = 1:height(resp)
            idx = find(elcNames == respLabels(i));
            if numel(idx) > 1
                error('Multiple electrodes rows match label "%s" (resp row %d).', respLabels(i), i);
            elseif numel(idx) == 1
                for j = 1:numel(validVarsToCopy)
                    vn = validVarsToCopy{j};
                    resp.(vn)(i) = electrodes.(vn)(idx);
                end
            end
        end














    end

    %% initalize response computation

    % Initialize times_unwarped and times columns in trials table
    trials.times_unwarped = cell(ntrials, 1);
    trials.times = cell(ntrials, 1);
    % Initialize columns for named epochs in response table (excluding 'base' since it populates resp.base)
    epochvars = epochs.Properties.RowNames(~ismember(epochs.Properties.RowNames, {'base'}));
    for iepoch = 1:length(epochvars)
        thisepoch = epochvars{iepoch};
        resp{:, thisepoch} = cel_chans_trials_nan;
    end
    % Ensure 'early_overlap_allowed' column exists (default to false if not provided)
    if ~ismember('early_overlap_allowed', epochs.Properties.VariableNames)
        epochs.early_overlap_allowed = false(nepochs, 1);
    end
    % Initialize X_epoch_early tracking columns in trials table ONLY where early overlap is allowed
    for iep = 1:nepochs
        ep_name = epochs.epoch{iep};
        if ~strcmp(ep_name, 'base') && epochs.early_overlap_allowed(iep)
            trials.([ep_name '_epoch_early']) = zeros(ntrials, 1);
        end
    end

    %% Pre-calculate unusable trials mask & combine with trials_to_analyze
    if ismember('unusable_trial', trials.Properties.VariableNames)
        if iscell(trials.unusable_trial)
            is_unusable_trial = cellfun(@(x) ~isempty(x) && logical(x), trials.unusable_trial);
        else
            unuse = trials.unusable_trial;
            unuse(isnan(unuse)) = 0; % assume nan indicates trial is not unusable
            is_unusable_trial = logical(unuse);
        end
    else
        is_unusable_trial = false(ntrials, 1);
    end

    % Trials to skip are either marked unusable or excluded by op.trials_to_analyze
    skip_trial = is_unusable_trial | ~op.trials_to_analyze;

%% 1. Pre-calculate all actual onset and offset times for each trial and check sequence order
    onset_times = zeros(ntrials, nepochs);
    offset_times = zeros(ntrials, nepochs);
    
    for itrial = 1:ntrials
        if skip_trial(itrial)
            continue;
        end
        for iepoch = 1:nepochs
            onset_times(itrial, iepoch) = parse_epoch_time(epochs.onset{iepoch}, trials, itrial);
            offset_times(itrial, iepoch) = parse_epoch_time(epochs.offset{iepoch}, trials, itrial);
        end

        % Skip trials with missing/NaN event timestamps
        if any(isnan(onset_times(itrial, :))) || any(isnan(offset_times(itrial, :)))
            skip_trial(itrial) = true;
            continue;
        end
        
        % Check that epochs appear in the order specified in the epochs table
        for iepoch = 1:(nepochs - 1)
            if onset_times(itrial, iepoch) > onset_times(itrial, iepoch + 1)
                if ~epochs.early_overlap_allowed(iepoch + 1)
                    error('Epoch order violation on trial %d: "%s" (onset %.4f) occurs after "%s" (onset %.4f), violating the specified epochs table order.', ...
                        itrial, epochs.epoch{iepoch}, onset_times(itrial, iepoch), ...
                        epochs.epoch{iepoch+1}, onset_times(itrial, iepoch+1));
                end
            end
        end
    end

    %% 2. Resolve Overlaps, Gaps, and Early Responses using Canonical Order
    % Omit NaNs to compute proper average chronological order across valid trials
    mean_onsets_all = mean(onset_times(~skip_trial, :), 1, 'omitnan');
    [~, chrono_order] = sort(mean_onsets_all);
    
    dt = mean(diff(D_in.time{1}));
    for itrial = 1:ntrials
        if skip_trial(itrial)
            continue;
        end
        res_onsets = onset_times(itrial, chrono_order);
        res_offsets = offset_times(itrial, chrono_order);
        allowed = epochs.early_overlap_allowed(chrono_order);
        names = epochs.epoch(chrono_order);
        
        for k = 1:nepochs
            % Clamp inherent negative durations
            if res_offsets(k) < res_onsets(k)
                res_offsets(k) = res_onsets(k);
            end
            
            % Forward-pass overlap check
            for j = (k + 1):nepochs
                if res_onsets(j) < res_offsets(k) - 1e-5
                    if allowed(j)
                        overlap_dur = res_offsets(k) - res_onsets(j);
                        trials{itrial, [names{j} '_epoch_early']} = ...
                            trials{itrial, [names{j} '_epoch_early']} + overlap_dur;
                        res_offsets(k) = res_onsets(j);
                        if res_onsets(j) < res_onsets(k)
                            res_onsets(k) = res_onsets(j);
                        end
                        warning('Overlap Allowed: Trial %d epoch "%s" started %.4fs early. Truncating preceding "%s" epoch.', ...
                            itrial, names{j}, overlap_dur, names{k});
                    else
                        error('Overlap detected on trial %d between epoch "%s" (ends at %.4f) and epoch "%s" (starts at %.4f).', ...
                            itrial, names{k}, res_offsets(k), names{j}, res_onsets(j));
                    end
                end
            end
        end
        % Check for gaps between adjacent epochs in canonical order exceeding average timestep duration
        for k = 1:(nepochs - 1)
            gap = res_onsets(k+1) - res_offsets(k);
            if gap > dt
                error('Gap detected on trial %d between epoch "%s" (ends at %.4f) and epoch "%s" (starts at %.4f), exceeding timestep duration (%.4fs).', ...
                    itrial, names{k}, res_offsets(k), names{k+1}, res_onsets(k+1), dt);
            end
        end
        
        onset_times(itrial, chrono_order) = res_onsets;
        offset_times(itrial, chrono_order) = res_offsets;
    end

    %% 3. Extract unwarped epoch-related responses and trial timecourses
    chrono_epochs = epochs(chrono_order, :);
    
    % Pre-compute 'base' responses first so they are available for baselining any epoch
    for itrial = 1:ntrials
        if skip_trial(itrial)
            continue;
        end
        for iepoch = 1:nepochs
            if strcmp(epochs.epoch{iepoch}, 'base')
                t_on = onset_times(itrial, iepoch);
                t_off = offset_times(itrial, iepoch);
                match_time_inds = D_in.time{1} > t_on & D_in.time{1} < t_off; 
                for ichan = 1:nchans
                    resp.base{ichan}(itrial) = mean(D_in.trial{1}(ichan, match_time_inds), 'includenan');
                end
            end
        end
    end
    for itrial = 1:ntrials
        if skip_trial(itrial)
            continue;
        end
        t_first = onset_times(itrial, chrono_order(1));
        t_last = offset_times(itrial, chrono_order(end));
        match_tr_inds = D_in.time{1} >= t_first & D_in.time{1} <= t_last;
        trials.times_unwarped{itrial} = D_in.time{1}(match_tr_inds);
        for iepoch = 1:nepochs
            thisep = epochs.epoch{iepoch};
            t_on = onset_times(itrial, iepoch);
            t_off = offset_times(itrial, iepoch);
            match_time_inds = D_in.time{1} > t_on & D_in.time{1} < t_off; 
            for ichan = 1:nchans
                if ~strcmp(thisep, 'base')
                    if iepoch == chrono_order(1)
                        cfg = [];
                        cfg.baseval = resp.base{ichan}(itrial); 
                        cfg.method = op.baseline_method; 
                        
                        tc_full = do_baselining(D_in.trial{1}(ichan, match_tr_inds), cfg); 
                        resp.timecourse_unwarped{ichan}{itrial} = tc_full;
                        if isnan(resp.base{ichan}(itrial)) || isempty(tc_full) || max(tc_full) > op.max_timecourse_base_ratio
                            resp.timecourse_unwarped{ichan}{itrial} = nan(size(tc_full));
                            resp.good_trial{ichan}(itrial) = false;
                        else
                            resp.good_trial{ichan}(itrial) = true;
                        end
                    end
                    if resp.good_trial{ichan}(itrial)
                        cfg = [];
                        cfg.baseval = resp.base{ichan}(itrial); 
                        cfg.method = op.baseline_method; 
                        
                        resp{ichan, thisep}{1}(itrial) = mean(do_baselining(D_in.trial{1}(ichan, match_time_inds), cfg)); 
                    else
                        resp{ichan, thisep}{1}(itrial) = nan; 
                    end
                end
            end
        end
    end

    %% 4. Calculate fixed segment sample counts for all epochs
    fs = 1 / mean(diff(D_in.time{1}));
    segment_dur_fix = chrono_epochs.dur_fix;
    segment_N = round(segment_dur_fix * fs);
    
    target_N_total = round(sum(segment_dur_fix) * fs);
    diff_N = target_N_total - sum(segment_N);
    if diff_N ~= 0
        [~, max_fill_idx] = max(segment_dur_fix);
        segment_N(max_fill_idx) = segment_N(max_fill_idx) + diff_N;
    end

    %% 5. Linearly timewarp timecourses back-to-back and build warped times vector
    for itrial = 1:ntrials
        if skip_trial(itrial)
            for ichan = 1:nchans
                resp.timecourse{ichan}{itrial} = nan(1, target_N_total);
            end
            continue;
        end
        t_first = onset_times(itrial, chrono_order(1));
        
        warped_time_chunks = cell(1, nepochs);
        running_t = t_first;
        
        for k = 1:nepochs
            nPts = segment_N(k);
            if nPts <= 0
                warped_time_chunks{k} = [];
                continue;
            end
            dur = segment_dur_fix(k);
            t_end_w = running_t + dur;
            
            warped_time_chunks{k} = linspace(running_t, t_end_w, nPts);
            running_t = t_end_w;
        end
        trials.times{itrial} = [warped_time_chunks{:}];
        for ichan = 1:nchans
            tc_unwarped = resp.timecourse_unwarped{ichan}{itrial};
            t_unwarped = trials.times_unwarped{itrial};
            
            if isempty(tc_unwarped) || all(isnan(tc_unwarped))
                resp.timecourse{ichan}{itrial} = nan(1, target_N_total);
                continue;
            end
            
            warped_chunks = cell(1, nepochs);
            
            for k = 1:nepochs
                nPts = segment_N(k);
                if nPts <= 0
                    warped_chunks{k} = [];
                    continue;
                end
                
                t_start = onset_times(itrial, chrono_order(k));
                t_end = offset_times(itrial, chrono_order(k));
                
                seg_mask = t_unwarped >= t_start & t_unwarped <= t_end;
                t_seg_actual = t_unwarped(seg_mask);
                tc_seg_actual = tc_unwarped(seg_mask);
                
                if length(t_seg_actual) < 2 || all(isnan(tc_seg_actual))
                    warped_chunks{k} = nan(1, nPts);
                else
                    t_target = linspace(t_seg_actual(1), t_seg_actual(end), nPts);
                    warped_chunks{k} = interp1(t_seg_actual, tc_seg_actual, t_target, 'linear', 'extrap');
                end
            end
            
            resp.timecourse{ichan}{itrial} = [warped_chunks{:}];
        end
    end

    resp.bad_elc = cellfun(@(x) all(isnan(x)), resp.base);
    resp.n_good_trials = cellfun(@(x) nnz(cellfun(@(y) ~all(isnan(y)), x)), resp.timecourse);
    if ~op.keep_unwarped_timecourse
        resp.timecourse_unwarped = [];
    end
end

%% Helper functions
function t = parse_epoch_time(spec, trial_table, itrial)
    if iscell(spec)
        evt_name = spec{1};
        shift = spec{2};
    else
        evt_name = spec;
        shift = 0;
    end
    t = trial_table{itrial, evt_name} + shift;
end

function normed_response = do_baselining(response, cfg)
    switch cfg.method
        case 'subtract'
            normed_response = response - cfg.baseval; 
        case 'subtract_then_divide'
            normed_response = (response - cfg.baseval) / cfg.baseval; 
    end
end