function resp = crossnobis_by_condition(resp, trials, cfg)
% CROSSNOBIS_BY_CONDITION Computes cross-validated Mahalanobis distance 
% (Crossnobis) for target variables within condition subsets per electrode.
%
% Inputs:
%   resp   - Table containing per-electrode response data & timecourses
%   trials - Table containing trial-level metadata (conditions, targets)
%   cfg    - Struct containing configuration parameters:
%            .target_var                   (string) variable in trials to decode (e.g. 'word')
%            .trial_condition_for_sorting  (string) variable in trials to split analysis (e.g. 'learn_con')
%            .sort_cond_vals               (cellstr/array) specific values of sort condition to use and their order (default [])
%            .epochs                       (table) epochs definition table (optional)
%            .epochs_to_analyze            (cellstr) epoch names matching epochs row names
%            .output_vars                  (cellstr) desired outputs: {'mean_dist','se_dist','snr','p_val','rdm'}
%            .proportion_of_max_folds      (double) proportion of trial pairs to sample (0..1, default 1)
%            .n_perm                       (int) number of permutations (default 0)
%            .cov_method                   (string) 'shrinkage' (default) or 'full'
%            .min_trials_per_sort_cond     (int) min good trials for sort condition (default 10)
%            .min_trials_per_target_label  (int) min good trials per target label (default 2)
%
% Output:
%   resp   - Updated response table with added cell variable 'cnob'

%% 1. Input Validation & Defaults
if ~isfield(cfg, 'target_var') || isempty(cfg.target_var)
    error('cfg.target_var must be specified (e.g., ''word'').');
end
if ~ismember(cfg.target_var, trials.Properties.VariableNames)
    error('Target variable ''%s'' not found in trials table.', cfg.target_var);
end

if ~isfield(cfg, 'trial_condition_for_sorting') || isempty(cfg.trial_condition_for_sorting)
    error('cfg.trial_condition_for_sorting must be specified (e.g., ''learn_con'').');
end
if ~ismember(cfg.trial_condition_for_sorting, trials.Properties.VariableNames)
    error('Sorting condition variable ''%s'' not found in trials table.', cfg.trial_condition_for_sorting);
end

% Sorting Condition Values & Order Option
if ~isfield(cfg, 'sort_cond_vals')
    cfg.sort_cond_vals = [];
end

% Validate Epochs Configuration
has_epochs_tbl = isfield(cfg, 'epochs') && ~isempty(cfg.epochs) && istable(cfg.epochs);
has_epochs_to_analyze = isfield(cfg, 'epochs_to_analyze') && ~isempty(cfg.epochs_to_analyze);

if has_epochs_to_analyze && ~has_epochs_tbl
    error('cfg.epochs_to_analyze was specified, but cfg.epochs table was not provided or is invalid.');
end

if ~has_epochs_tbl
    fprintf('cfg.epochs not provided. Analyzing full trial timecourse for each trial.\n');
    epochs_list = {'full_trial'};
else
    if isempty(cfg.epochs.Properties.RowNames)
        if ismember('epoch', cfg.epochs.Properties.VariableNames)
            cfg.epochs.Properties.RowNames = cfg.epochs.epoch;
        else
            error('cfg.epochs table must have RowNames set or a variable named ''epoch''.');
        end
    end
    
    if has_epochs_to_analyze
        missing_epochs = setdiff(cfg.epochs_to_analyze, cfg.epochs.Properties.RowNames);
        if ~isempty(missing_epochs)
            error('The following epochs_to_analyze were not found in cfg.epochs: %s', strjoin(missing_epochs, ', '));
        end
        epochs_list = cfg.epochs_to_analyze(:)';
    else
        epochs_list = cfg.epochs.Properties.RowNames(:)';
    end
end

% Validate Output Variables
valid_outputs = {'mean_dist', 'se_dist', 'snr', 'p_val', 'rdm'};
if ~isfield(cfg, 'output_vars') || isempty(cfg.output_vars)
    cfg.output_vars = valid_outputs;
else
    invalid_vars = setdiff(cfg.output_vars, valid_outputs);
    if ~isempty(invalid_vars)
        error('Invalid output variable(s) in cfg.output_vars: %s. Valid choices: %s', ...
            strjoin(invalid_vars, ', '), strjoin(valid_outputs, ', '));
    end
end

% Trial Threshold Defaults
if ~isfield(cfg, 'min_trials_per_sort_cond') || isempty(cfg.min_trials_per_sort_cond)
    cfg.min_trials_per_sort_cond = 10;
end
if ~isfield(cfg, 'min_trials_per_target_label') || isempty(cfg.min_trials_per_target_label)
    cfg.min_trials_per_target_label = 2;
end

% General Parameters Defaults
if ~isfield(cfg, 'proportion_of_max_folds') || isempty(cfg.proportion_of_max_folds)
    cfg.proportion_of_max_folds = 1.0;
end
cfg.proportion_of_max_folds = max(0, min(1, cfg.proportion_of_max_folds));

if ~isfield(cfg, 'n_perm') || isempty(cfg.n_perm)
    cfg.n_perm = 0;
end

if ~isfield(cfg, 'cov_method') || isempty(cfg.cov_method)
    cfg.cov_method = 'shrinkage';
end

%% 2. Setup Metadata
nchans = height(resp);

% Determine condition values and their ordering
if ~isempty(cfg.sort_cond_vals)
    cond_vals = cfg.sort_cond_vals;
    if ischar(cond_vals)
        cond_vals = {cond_vals};
    end
else
    cond_vals = unique(trials.(cfg.trial_condition_for_sorting), 'stable');
    if iscell(cond_vals)
        cond_vals = cond_vals(~cellfun(@(x) any(ismissing(x)), cond_vals));
    else
        cond_vals = cond_vals(~isnan(cond_vals));
    end
end

nconds = length(cond_vals);
nepochs = length(epochs_list);

% Full set of unique target labels across all trials (defines fixed RDM size)
all_unique_targets = unique(trials.(cfg.target_var));
if iscell(all_unique_targets)
    all_unique_targets = all_unique_targets(~cellfun(@(x) any(ismissing(x)), all_unique_targets));
else
    all_unique_targets = all_unique_targets(~isnan(all_unique_targets));
end

cnob_cell = cell(nchans, 1);

%% 3. Main Processing Loop
for ichan = 1:nchans
    if ismember('good_trial', resp.Properties.VariableNames) && ~isempty(resp.good_trial{ichan})
        good_trials_chan = logical(resp.good_trial{ichan});
    else
        good_trials_chan = true(height(trials), 1);
    end
    
    if ismember('bad_elc', resp.Properties.VariableNames) && resp.bad_elc(ichan)
        cnob_cell{ichan} = [];
        continue;
    end

    epoch_results_table = table('Size', [nepochs, 1], ...
        'VariableTypes', {'cell'}, ...
        'VariableNames', {'trialconds'}, ...
        'RowNames', epochs_list);

    for iepoch = 1:nepochs
        this_epoch = epochs_list{iepoch};
        X_trials = extract_epoch_data(resp, ichan, this_epoch);
        
        if isempty(X_trials)
            continue;
        end

        cond_rows = cell(nconds, 1);
        
        for icond = 1:nconds
            if iscell(cond_vals) || isstring(cond_vals)
                this_cond = cond_vals{icond};
            else
                this_cond = cond_vals(icond);
            end

            if iscell(trials.(cfg.trial_condition_for_sorting)) || isstring(trials.(cfg.trial_condition_for_sorting))
                cond_mask = strcmp(trials.(cfg.trial_condition_for_sorting), this_cond) & good_trials_chan;
            else
                cond_mask = (trials.(cfg.trial_condition_for_sorting) == this_cond) & good_trials_chan;
            end
            
            y_sub = trials.(cfg.target_var)(cond_mask);
            X_sub = X_trials(cond_mask, :);
            
            % Remove missing/NaN data
            if iscell(y_sub)
                nan_trials = any(isnan(X_sub), 2) | cellfun(@(x) any(ismissing(x)), y_sub);
            else
                nan_trials = any(isnan(X_sub), 2) | isnan(y_sub);
            end
            X_sub(nan_trials, :) = [];
            y_sub(nan_trials) = [];

            % Check 1: Sort Condition Minimum Total Trials
            if size(X_sub, 1) < cfg.min_trials_per_sort_cond
                cond_struct = init_empty_results(cfg.output_vars, all_unique_targets);
            else
                % Check 2: Filter target labels meeting min_trials_per_target_label
                [valid_y_sub, valid_X_sub] = filter_labels_by_count(...
                    y_sub, X_sub, cfg.min_trials_per_target_label);
                
                unique_valid_labels = unique(valid_y_sub);
                
                if length(unique_valid_labels) < 2
                    cond_struct = init_empty_results(cfg.output_vars, all_unique_targets);
                else
                    % Compute Crossnobis distance
                    [mean_d, se_d, snr_val, rdm_mat] = compute_crossnobis(...
                        valid_X_sub, valid_y_sub, all_unique_targets, ...
                        cfg.proportion_of_max_folds, cfg.cov_method);

                    % Permutation testing
                    if cfg.n_perm > 0 && ismember('p_val', cfg.output_vars)
                        p_val_num = run_permutation_test(...
                            valid_X_sub, valid_y_sub, all_unique_targets, mean_d, ...
                            cfg.n_perm, cfg.proportion_of_max_folds, cfg.cov_method);
                    else
                        p_val_num = NaN;
                    end

                    % Assemble output struct (wrap rdm matrix in cell for 1-row table conversion)
                    cond_struct = struct();
                    if ismember('mean_dist', cfg.output_vars), cond_struct.mean_dist = mean_d; end
                    if ismember('se_dist',   cfg.output_vars), cond_struct.se_dist   = se_d;   end
                    if ismember('snr',       cfg.output_vars), cond_struct.snr       = snr_val; end
                    if ismember('p_val',     cfg.output_vars), cond_struct.p_val     = p_val_num; end
                    if ismember('rdm',       cfg.output_vars), cond_struct.rdm       = {rdm_mat}; end
                end
            end
            
            cond_rows{icond} = struct2table(cond_struct);
        end

        trialconds_tbl = vertcat(cond_rows{:});
        trialconds_tbl.Properties.RowNames = string(cond_vals);
        epoch_results_table.trialconds{iepoch} = trialconds_tbl;
    end

    cnob_cell{ichan} = epoch_results_table;
end

resp.cnob = cnob_cell;
end

%% Helper Functions

function [y_out, X_out] = filter_labels_by_count(y, X, min_count)
    unique_y = unique(y);
    keep_mask = false(size(y));
    for i = 1:length(unique_y)
        if iscell(y)
            idx = strcmp(y, unique_y{i});
        else
            idx = (y == unique_y{i});
        end
        if sum(idx) >= min_count
            keep_mask(idx) = true;
        end
    end
    y_out = y(keep_mask);
    X_out = X(keep_mask, :);
end

function X_trials = extract_epoch_data(resp, ichan, epoch_name)
    X_trials = [];
    if ismember(epoch_name, resp.Properties.VariableNames)
        raw_data = resp{ichan, epoch_name};
        if iscell(raw_data)
            if numel(raw_data) == 1
                inner_data = raw_data{1};
                if isnumeric(inner_data)
                    X_trials = inner_data;
                elseif iscell(inner_data)
                    X_trials = cell2mat(inner_data(:));
                end
            else
                if isnumeric(raw_data{1})
                    X_trials = cell2mat(raw_data(:));
                elseif iscell(raw_data{1})
                    X_trials = cell2mat(raw_data{1}(:));
                end
            end
        elseif isnumeric(raw_data)
            X_trials = raw_data;
        end
    elseif ismember('timecourse', resp.Properties.VariableNames)
        tc = resp.timecourse{ichan};
        if iscell(tc)
            if numel(tc) == 1 && iscell(tc{1})
                X_trials = cell2mat(tc{1}(:));
            else
                X_trials = cell2mat(tc(:));
            end
        elseif isnumeric(tc)
            X_trials = tc;
        end
    end
end

function [mean_d, se_d, snr_val, rdm_mat] = compute_crossnobis(X, y, all_targets, prop_folds, cov_method)
    K_all = length(all_targets);
    rdm_mat = nan(K_all, K_all);
    
    unique_present = unique(y);
    K_present = length(unique_present);
    
    % Residual noise covariance calculation
    X_res = zeros(size(X));
    for k = 1:K_present
        if iscell(y)
            k_idx = strcmp(y, unique_present{k});
        else
            k_idx = (y == unique_present{k});
        end
        X_res(k_idx, :) = X(k_idx, :) - mean(X(k_idx, :), 1);
    end
    
    if strcmp(cov_method, 'shrinkage')
        Sigma = cov_shrinkage(X_res);
    else
        Sigma = cov(X_res);
    end
    Sigma_inv = pinv(Sigma);

    pair_dists_all = [];

    for a = 1:K_all
        target_a = all_targets{a};
        if iscell(y)
            idx_a = find(strcmp(y, target_a));
        else
            idx_a = find(y == target_a);
        end
        if length(idx_a) < 2, continue; end
        
        for b = (a+1):K_all
            target_b = all_targets{b};
            if iscell(y)
                idx_b = find(strcmp(y, target_b));
            else
                idx_b = find(y == target_b);
            end
            if length(idx_b) < 2, continue; end
            
            Na = length(idx_a);
            Nb = length(idx_b);
            
            [grid_a, grid_b] = meshgrid(1:Na, 1:Nb);
            all_pairs = [grid_a(:), grid_b(:)];
            n_total_pairs = size(all_pairs, 1);
            
            n_sample_pairs = max(1, ceil(prop_folds * n_total_pairs));
            if n_sample_pairs < n_total_pairs
                sample_idx = randperm(n_total_pairs, n_sample_pairs);
                pairs_to_run = all_pairs(sample_idx, :);
            else
                pairs_to_run = all_pairs;
            end
            
            pair_dists = zeros(n_sample_pairs, 1);
            
            for p = 1:n_sample_pairs
                ia = idx_a(pairs_to_run(p, 1));
                ib = idx_b(pairs_to_run(p, 2));
                
                % Test difference vector
                d_test = X(ia, :) - X(ib, :);
                
                % Training difference vector (excluding test pair)
                mu_a_train = mean(X(setdiff(idx_a, ia), :), 1);
                mu_b_train = mean(X(setdiff(idx_b, ib), :), 1);
                d_train = mu_a_train - mu_b_train;
                
                % Crossnobis dot product
                pair_dists(p) = d_train * Sigma_inv * d_test';
            end
            
            d_ab = mean(pair_dists);
            rdm_mat(a, b) = d_ab;
            rdm_mat(b, a) = d_ab;
            pair_dists_all = [pair_dists_all; pair_dists];
        end
    end
    
    valid_dists = rdm_mat(triu(true(K_all), 1) & ~isnan(rdm_mat));
    if isempty(valid_dists)
        mean_d = NaN;
        se_d = NaN;
        snr_val = NaN;
    else
        mean_d = mean(valid_dists);
        se_d = std(pair_dists_all) / sqrt(length(pair_dists_all));
        if se_d > 0
            snr_val = mean_d / se_d;
        else
            snr_val = 0;
        end
    end
end

function p_val = run_permutation_test(X, y, all_targets, actual_mean_d, n_perm, prop_folds, cov_method)
    if isnan(actual_mean_d)
        p_val = NaN;
        return;
    end
    null_means = zeros(n_perm, 1);
    for ip = 1:n_perm
        y_shuffled = y(randperm(length(y)));
        [null_mean_d, ~, ~, ~] = compute_crossnobis(X, y_shuffled, all_targets, prop_folds, cov_method);
        null_means(ip) = null_mean_d;
    end
    p_val = (sum(null_means >= actual_mean_d) + 1) / (n_perm + 1);
end

function Sigma = cov_shrinkage(X)
    [N, P] = size(X);
    S = (X' * X) / N;
    prior = trace(S) / P * eye(P);
    
    d2 = norm(S - prior, 'fro')^2;
    if d2 == 0
        Sigma = S;
        return;
    end
    
    y = zeros(P, P, N);
    for i = 1:N
        y(:, :, i) = X(i, :)' * X(i, :);
    end
    
    b2 = 0;
    for i = 1:N
        b2 = b2 + norm(y(:, :, i) - S, 'fro')^2;
    end
    b2 = b2 / (N^2);
    
    shrinkage = max(0, min(1, b2 / d2));
    Sigma = (1 - shrinkage) * S + shrinkage * prior;
end

function empty_struct = init_empty_results(output_vars, all_unique_targets)
    empty_struct = struct();
    K = length(all_unique_targets);
    if ismember('mean_dist', output_vars), empty_struct.mean_dist = NaN; end
    if ismember('se_dist',   output_vars), empty_struct.se_dist   = NaN; end
    if ismember('snr',       output_vars), empty_struct.snr       = NaN; end
    if ismember('p_val',     output_vars), empty_struct.p_val     = NaN; end
    if ismember('rdm',       output_vars), empty_struct.rdm       = {nan(K, K)}; end
end