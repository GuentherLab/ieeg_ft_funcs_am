function [resp, trials] = get_epoched_responses(D_in,trials,op)


    epochs = op.epochs; 
    nepochs = height(epochs); 

    % table containing responses during epochs for each chan
    ntrials = height(trials);
    nchans = length(D_in.label);
    nans_tr = nan(ntrials,1); 
    false_tr = false(ntrials,1);
    cel_trials = cell(ntrials,1); 
    cel_chans_trials_nan = repmat({nans_tr},nchans,1); % 1 element per trial per chan
    cel_chans_trials_false = repmat({false_tr},nchans,1);
    resp = table(    D_in.label, cel_chans_trials_nan,   repmat({cel_trials},nchans,1), cel_chans_trials_false, ....
      'VariableNames', {'chan', 'base',                  'timecourse',                  'good_trial'            }); 
    epochvars = epochs.Properties.RowNames(~ismember(epochs.Properties.RowNames, {'trial', 'base'}));
    for iepoch = 1:length(epochvars)
        thisepoch = epochvars{iepoch};
        resp{:,thisepoch} = cel_chans_trials_nan;
    end


    % make sure 'base' and 'trial' are on top; base required first so we can do baselining on all other epochs
    assert([ismember('trial',epochs.Properties.RowNames) & ismember('base',epochs.Properties.RowNames)], ['''trial'' and ''base'' must be rows in op.epochs'])
    [~, idx] = ismember(epochs.Properties.RowNames, {'base', 'trial'});
            epochs = [epochs(idx == 1, :); epochs(idx == 2, :); epochs(idx == 0, :)]; 
    

    
    % extract epoch-related responses, get phonemes on each trial
    %%%% trials.times{itrial} use global time coordinates
    %%%% ....... start at a fixed baseline window before stim onset
    %%%% ....... end at a fixed time buffer after speech offset
    for itrial = 1:ntrials % itrial is absolute index across sessions; does not equal "trial_id" from loaded tables
        for iepoch = 1:nepochs
            thisep = epochs.epoch{iepoch};

            % process epoch onset and offset specifications
            %%%% these can either be an event name in the trial trial, or a name and timeshift
            thisep_onset = epochs.onset{iepoch};
            if length(thisep_onset) == 1 % if only an event name was specified
                onset_name = thisep_onset;
                onset_shift = 0; % if no shifting off of the named event was specified, set it to zero
            elseif length(thisep_onset) == 2 % if both an event name and time shift was specified
                onset_name = thisep_onset{1};
                onset_shift = thisep_onset{2};
            end

            thisep_offset = epochs.offset{iepoch};
            if length(thisep_offset) == 1 % if only an event name was specified
                offset_name = thisep_offset;
                offset_shift = 0; % if no shifting off of the named event was specified, set it to zero
            elseif length(thisep_offset) == 2  % if both an event name and time shift was specified
                offset_name = thisep_offset{1};
                offset_shift = thisep_offset{2};
            end
            assert(all(ismember({onset_name; offset_name}, trials.Properties.VariableNames)), ['epoch event name not found in trial table'])
            assert(all(isnumeric([onset_shift; offset_shift])),'epoch event shift not numeric')


            % find time inds constituting this epoch in this trial
            match_time_inds = D_in.time{1} > trials{itrial,onset_name}+onset_shift & D_in.time{1} < trials{itrial,offset_name}+offset_shift; 

            % if epoch is named 'trial', add 'starts' and 'ends' variables to trialtable and use it to create 'timecourse' variable, but don't process it further
            % if epoch is named 'base', use it to compute baselines for other epochs
            if strcmp(thisep,'trial')
                trials.times{itrial} = D_in.time{1}(match_time_inds); % times in this trial window
            end

            % get response for each chan in each epoch
            for ichan = 1:nchans
                switch thisep
                    case 'base'
                        % use mean rather than nanmean, so that trials which had artifacts marked with NaNs will be excluded
                        resp.base{ichan}(itrial) = mean( D_in.trial{1}(ichan, match_time_inds), 'includenan' ); % mean response during baseline

                    case 'trial'
                        cfg = [];
                        cfg.baseval = resp.base{ichan}(itrial); 
                        cfg.method = op.baseline_method; 
                        resp.timecourse{ichan}{itrial} = do_baselining(D_in.trial{1}(ichan, match_time_inds), cfg); 

                       %%% if response looks artifactually high or if baseline is nan, set/leave all response values for this trial to nan
                       if isnan(resp.base{ichan}(itrial))   ||   max(resp.timecourse{ichan}{itrial}) > op.max_timecourse_base_ratio
                            resp.timecourse{ichan}{itrial} = nan(size(resp.timecourse{ichan}{itrial}));
                             resp.good_trial{ichan}(itrial) = false;
                       else
                           resp.good_trial{ichan}(itrial) = true;
                       end

                    otherwise %%% all other epochs
                        if resp.good_trial{ichan}(itrial)
                            cfg = [];
                            cfg.baseval = resp.base{ichan}(itrial); 
                            cfg.method = op.baseline_method; 
                            resp{ichan,thisep}{1}(itrial) = mean(do_baselining(D_in.trial{1}(ichan, match_time_inds), cfg)); 
                        elseif ~resp.good_trial{ichan}(itrial)
                            resp{ichan,thisep}{1}(itrial) = nan; 
                        end
                end
            end
        end
    end
    resp.bad_elc = cellfun(@(x)all(isnan(x)),resp.base);
end



%% subfunction - takes a response (numerical array) and does baseline normalization used a specified method
function normed_response = do_baselining(response,cfg)
    switch cfg.method
        case 'subtract'
            normed_response = response - cfg.baseval; 

        case 'subtract_then_divide'
            normed_response = [response - cfg.baseval] / cfg.baseval; 
    end
end

