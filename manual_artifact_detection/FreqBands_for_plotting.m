% window:
triplet_subjects = [];
SMSL_subjects = ["sub-DM1005", "sub-DM1007", "sub-DM1008", "sub-DM1024", "sub-DM1025", "sub-DM1037"];

%% triplet loop
% store sessions as rows, and trials as columns

%for sub=1:length(triplet_subjects)
    %SUBJECT = triplet_subjects(i); 
    SUBJECT = 'DBS3001';
    disp(SUBJECT);

    clear continuous_FT
    clear D_wavtransf
    clear trial_startend
    clear timepoints
    clear trialed

    PATH_DATASET = 'Z:\DBS';
    PATH_SUB = [PATH_DATASET filesep SUBJECT];
    PATH_PREPROC_DATA = [PATH_SUB filesep 'Preprocessed Data'];
    PATH_FT = [PATH_PREPROC_DATA filesep 'FieldTrip'];
    PATH_SYNC = [PATH_PREPROC_DATA filesep 'Sync'];
    PATH_ANNOT = [PATH_SYNC filesep 'annot'];

    filename = [SUBJECT '_ft_raw_session.mat'];

    tempname = load([PATH_FT filesep filename]);
    FT_file = tempname.D;

    % create new table with .trial and .time stored continuously
    %continuous_FT = struct('label','trial','time','cfg','hdr');
    %{
    continuous_FT.label = FT_file.label;
    continuous_FT.cfg = FT_file.cfg;
    continuous_FT.hdr = FT_file.hdr;
    continuous_FT.trial = {};
    continuous_FT.trial{1,1} = horzcat(FT_file.trial{1,1},...
                                    FT_file.trial{1,2},...
                                    FT_file.trial{1,3},...
                                    FT_file.trial{1,4});
    continuous_FT.time = {};
    continuous_FT.time{1,1} = horzcat(FT_file.time{1,1},...
                                    FT_file.time{1,2},...
                                    FT_file.time{1,3},...
                                    FT_file.time{1,4});
    %}

    % loop through each of the sessions
    for j = 1:length(FT_file.trial)
        % size to trial matrix to determine size of transformed trial structure
        sz = size(FT_file.trial{1,j});
    
        fsample = round(1/mean(diff(FT_file.time{1})));
    
        envelope_wavpow_outfreq = 100;
        raw2wavpow_confactor = fsample/envelope_wavpow_outfreq;
        %new_sz = round(sz(2)/raw2wavpow_confactor);
        new_sz = floor(sz(2)/raw2wavpow_confactor);
        
        % calculate frequencies
        freqs = 10.^(linspace(0,2.3,20)); % number of frequencies: 20 (logarithmaclly spaced between 1 and 200 [previously between 1 and 316])
        D_wavtransf.trial{j,1} = zeros([sz(1),new_sz,length(freqs)]);
        for ifreq = 1:length(freqs)
            clear temp_in
            fprintf('ifreq: %d\n', ifreq);
    
            thisfreq = freqs(ifreq);
            cfg=[];
            cfg.out_freq = 100;
            cfg.wav_freq = thisfreq;
            cfg.wav_width = 7;

            % create a temporary ft file for bml_envelope_wavpow to use
            temp_in.trial{1,1} = FT_file.trial{1,j};
            temp_in.time{1,1} = FT_file.time{1,j};
            temp_in.fsample = fsample;
    
            %D_wavtransf{ifreq} = bml_envelope_wavpow(cfg,FT_file);
            temp_out = bml_envelope_wavpow(cfg,temp_in);
    
            D_wavtransf.trial{j,1}(:,:,ifreq) = temp_out.trial{1,1};
            D_wavtransf.time{j,1} = temp_out.time{1,1};
            D_wavtransf.cfg = temp_out.cfg;
            D_wavtransf.cfg = rmfield(D_wavtransf.cfg,'wav_freq');
        end
    end

    D_wavtransf.label = FT_file.label;
    D_wavtransf.hdr = FT_file.hdr;

    % save untrialed data
    %save([PATH_FT filesep SUBJECT '_ft_raw_session_freqbands.mat'],'D_wavtransf','-v7.3');

    % break up into trials
    size_D = size(D_wavtransf.trial{1,1});
    numFreq = size_D(3);
    annot_syllable = readtable([PATH_ANNOT filesep SUBJECT '_produced_syllable'],'Delimiter','\t');

    % create a list of only the end of the last syllable
    annot_end = annot_syllable.ends(3:3:end);

    % where should I get the trial window from for triplet subjects?
    annot_stim = bml_annot_read_tsv([PATH_ANNOT filesep SUBJECT '_stimulus_triplet']);

    disp('determining where the trials start and end');
    % loop through each of the sessions
    % row indicates trial
    % column indicates on or offset
    % z axis indicates session
    %trial_startend = zeros(length(annot_end),2,length(FT_file.trial));
    j = 1; % trial number
    prev_j = 1;
    for i=1:length(FT_file.trial) % loop through the sessions
        % determine where the trials start and end

        while_count = 1; % index for trial_startend
        %trial_startend(:,2) = annot.sp_off(:) + 1.5;
        while j <= length(annot_stim.session_id) && annot_stim.session_id(j) == i
            %j=j+1;
            rownum_thirdsyl = j*3; % row number for the third syllable in _produced_syllable
    
            if ~isnan(annot_end(j)) % third syllable offset
                % if the speech off of the third syllable exists, use that
                trial_startend(while_count,2,i) = annot_end(j) + 1.5;
                %disp('1');
            elseif ~isnan(annot_syllable.starts(rownum_thirdsyl)) % third syllable onset
                trial_startend(while_count,2,i) = annot_syllable.starts(rownum_thirdsyl) + 1.5;
                %disp('2');
            elseif ~isnan(annot_syllable.ends(rownum_thirdsyl-1)) % second syllable offset
                trial_startend(while_count,2,i) = annot_syllable.ends(rownum_thirdsyl-1) + 1.5;
                %disp('3');
            elseif ~isnan(annot_syllable.starts(rownum_thirdsyl-1)) % second syllable onset
                trial_startend(while_count,2,i) = annot_syllable.starts(rownum_thirdsyl-1) + 1.5;
                %disp('4');
            elseif ~isnan(annot_syllable.ends(rownum_thirdsyl-2)) % first syllable offset
                trial_startend(while_count,2,i) = annot_syllable.ends(rownum_thirdsyl-2) + 1.5;
                %disp('5');
            elseif ~ isnan(annot_syllable.starts(rownum_thirdsyl-2)) % first syllable onset
                trial_startend(while_count,2,i) = annot_syllable.starts(rownum_thirdsyl-2) + 1.5;
                %disp('6');
            else
                fprintf('no triplet speech onset and offset data, used stim offset. trial number: %d \n', j);
                trial_startend(while_count,2,i) = annot_stim.ends(j) + 4;
            end

            j=j+1;
            while_count = while_count + 1;
        end

        trial_startend(:,1,i) = annot_stim.starts(prev_j:(j-1)) - 0.5;
        prev_j = j;
        %fprintf('while break row num: %d\n', j);
    end
    
    disp('calculating timepoints');
    % find all timepoints between start and end of each trial
    timepoints = zeros(length(trial_startend),2,length(FT_file.trial)); % each row is a trial; first column is start location; second column is end location; z axis is session
        % both timepoints are from continuous
    for i=1:length(FT_file.trial) % looping through the sessions
        j=1;
        annot_sz = size(trial_startend(:,:,i));
        for k = 1:annot_sz(1) % looping through the trials
            % for each timepoint (a) iterate through continuous and for each column in continuous (b) do a-b
            % when a-b becomes < 0, the previous timepoint is the selected column
            % optimization: after one timepoint is calculated, the calculation does not need to be performed on the previous set of times in continuous
        
            % loop to determine onset column
            onset=1;
            if trial_startend(k,1,i) == 0 || trial_startend(k,2,i) == 0
                % if either equal 0, there is no data so the timepoints should be 0 as well
                timepoints(k,1,i) = 0;
                timepoints(k,2,i) = 0;
            else
                while onset<trial_startend(k,1,i)
                    onset = D_wavtransf.time{i,1}(1,j,1);
                    j=j+1;
                end
                timepoints(k,1,i) = j;
                
                stop=1;
                while stop<trial_startend(k,2,i)
                    stop = D_wavtransf.time{i,1}(1,j,1);
                    j=j+1;
                end
                timepoints(k,2,i) = j;
    
                j = 1; % timepoints may overlap, have to reset j
            end
        end
    end
    
    disp('breaking up into trials');
    %trialed = struct('label','trial','time','cfg','hdr','fsample');
    trialed.label = D_wavtransf.label;
    trialed.trial = {};
    trialed.time = {};
    trialed.cfg = D_wavtransf.cfg;
    trialed.hdr = D_wavtransf.hdr;
    %trialed.fsample = D_wavtransf.fsample;

    %trialed = cell(1,numFreq);
    for i=1:length(FT_file.trial) % looping through the sessions
        % copy other data stored in the table to the new variable (basically not .time and .trial) 
        
        % paste the data into the new table under the correct trial number
        % run through each trial

        % delete rows containing 0 (to account for differently sized sessions)
        temp = trial_startend(:,:,i);
        delRows = find(trial_startend(:,1,i) == 0 | trial_startend(:,2,i) == 0);
        temp(delRows,:) = [];
        
        num_trials = size(temp);
        for j=1:num_trials(1) % trial
            %trial_temp = [];
            %count = 1;
            %{
            for k=timepoints(j,1):timepoints(j,2) % timepoints
                % take each datapoint between onset and offset from continous and put into .trial
                % for each electrode

                %trial_sz = size(continuous_FT{1,n}.trial);
                trial_sz = size(D_wavtransf.trial{:,:,1});
                for m=1:trial_sz(1);
                    trial_temp(m,count) = continuous_FT{1,n}.trial(m,k); % need to do for each electrode
                end
                
                % take each timepoint between onset and offset and put into .time
                time_temp(1,count) = continuous_FT{1,n}.trial(1,k);
                count = count+1;
            end
            %}
            trialed.trial{i,j} = D_wavtransf.trial{i,1}(:,timepoints(j,1,i):timepoints(j,2,i),:);
            trialed.time{i,j} = D_wavtransf.time{i,1}(:,timepoints(j,1,i):timepoints(j,2,i));
        end     
    end

    save([PATH_FT filesep SUBJECT '_ft_raw_session_freqbands-trial.mat'],'trialed','-v7.3');
    disp('finished');
%end

%% SMSL loop
%for sub=1:length(SMSL_subjects)
    %SUBJECT = char(SMSL_subjects(i)); 
    SUBJECT = 'sub-DM1005';
    disp(SUBJECT);

    clear D_wavtransf
    clear trial_startend
    clear timepoints
    clear trialed

    PATH_DATASET = 'Y:\DBS';
    PATH_DER = [PATH_DATASET filesep 'derivatives'];
    PATH_SUB = [PATH_DER filesep SUBJECT];
    PATH_FT = [PATH_SUB filesep 'fieldtrip'];
    PATH_ANNOT = [PATH_SUB filesep 'annot'];

    filename = [SUBJECT '_ses-intraop_task-smsl_ft-raw'];

    tempname = load([PATH_FT filesep filename]);
    FT_file = tempname.D;

    % computes trial start and end time
    disp('determining where the trials start and end');
    %annot = readtable([PATH_ANNOT filesep SUBJECT '_ses-intraop_task-smsl_annot-trials.tsv'],'FileType','text','Delimiter','\t');
    annot = readtable([PATH_ANNOT filesep SUBJECT '_ses-intraop_task-smsl_annot-produced-syllables.tsv'],'FileType','text','Delimiter','\t');

    trial_startend(:,1) = annot.audio_onset(:) - 0.5;
    %trial_startend(:,2) = annot.sp_off(:) + 1.5;
    for j=1:length(annot.sp_off)
        if isnan(annot.sp_off(j)) && ~isnan(annot.keypress_time(j))
            trial_startend(j,2) = annot.keypress_time(j); % if the trial is unusable and there is a keypress time it will instead get data ending at the keypress time
        elseif isnan(annot.sp_off(j)) && isnan(annot.keypress_time(j))
            trial_startend(j,2) = (annot.audio_onset(j) - 0.5) + 5;
        else
            trial_startend(j,2) = annot.sp_off(j) + 1.5;
        end
    end

    % compute frequencies
    disp('computing frequencies');
    D_wavtransf.label = FT_file.label;
    D_wavtransf.hdr = FT_file.hdr;
    D_wavtransf.fsample = FT_file.fsample;
    D_wavtransf.trial = {};

    % size of trial matrix to determine size of transformed trial structure
    sz = size(FT_file.trial{1,1});

    envelope_wavpow_outfreq = 100;

    raw2wavpow_confactor = FT_file.fsample/envelope_wavpow_outfreq;
    new_sz = round(sz(2)/raw2wavpow_confactor);

    freqs = 10.^(linspace(0,2.3,20)); % number of frequencies: 20 (logarithmaclly spaced between 1 and 200 [previously between 1 and 316])
    D_wavtransf.trial{1,1} = zeros([sz(1),new_sz,length(freqs)]);
    for ifreq = 1:length(freqs)
        fprintf('ifreq: %d \n', ifreq);
        thisfreq = freqs(ifreq);
        cfg=[];
        cfg.out_freq = envelope_wavpow_outfreq;
        cfg.wav_freq = thisfreq;
        cfg.wav_width = 7;

        temp_out = bml_envelope_wavpow(cfg,FT_file);
        %D_wavtransf{ifreq} = bml_envelope_wavpow(cfg,FT_file); 

        D_wavtransf.trial{1,1}(:,:,ifreq) = temp_out.trial{1,1};
        D_wavtransf.time = temp_out.time;
        %{
        for iD = 1:sz(1)
            for jD = 1:sz(2)
                D_wavtransf.trial{1,1}(iD,jD,ifreq) = temp.trial{1,1}(iD,jD);
            end
        end
        %}
        D_wavtransf.cfg = temp_out.cfg;
        D_wavtransf.cfg = rmfield(D_wavtransf.cfg,'wav_freq');
    end

    % break up into trials
    numFreq = length(freqs);

    % do I need the continuous variable anymore?
    %{
    continuous = cell(1,numFreq);
    for j = 1:numFreq
        continuous{1,j} = struct('time','trial');
        continuous{1,j}.time = D_wavtransf{1,j}.time{:,:};
        continuous{1,j}.trial = D_wavtransf{1,j}.trial{:,:};
    end
    %}

    save([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_ft-raw_freqbands.mat'],'D_wavtransf','-v7.3');

    disp('calculating timepoints');
    % find all timepoints between start and end of each trial
    timepoints = []; % each row is a trial; first column is start location; second column is end location, locations column number
        % both timepoints are from continuous
    annot_sz = size(annot);
    j=1;
    for k = 1:annot_sz(1)
        % for each timepoint (a) iterate through continuous and for each column in continuous (b) do a-b
        % when a-b becomes < 0, the previous timepoint is the selected column
        % optimization: after one timepoint is calculated, the calculation does not need to be performed on the previous set of times in continuous
    
        % loop to determine onset column
        onset=1;
        while onset<trial_startend(k,1)
            %onset = continuous{1,1}.time(1,j);
            onset = D_wavtransf.time{1,1}(1,j);
            j=j+1;
        end
        timepoints(k,1) = j;
        
        stop=1;
        %offsettimep = annot{k,1} + annot{k,2}; % annot{i,2} is the duration, not offset
        while stop<trial_startend(k,2)
            %stop = continuous{1,1}.time(1,j);
            stop = D_wavtransf.time{1,1}(1,j);
            j=j+1;
        end
        timepoints(k,2) = j;

        j = 1; % the trials may overlap, so it is necessary to reset j
    end
    
    disp('breaking up into trials');
    
    %trialed = cell(1,numFreq);
    trialed = struct('label','trial','time','cfg','hdr','fsample');
    trialed.label = D_wavtransf.label;
    trialed.trial = {};
    trialed.time = {};
    trialed.cfg = D_wavtransf.cfg;
    trialed.hdr = D_wavtransf.hdr;
    trialed.fsample = D_wavtransf.fsample;
    % loop through each frequency
    
        % copy other data stored in the table to the new variable (basically not .time and .trial)
        %{
        trialed{1,n} = struct('label','trial','time','cfg','hdr','fsample');
        trialed{1,n}.label = D_wavtransf{1,n}.label;
        trialed{1,n}.trial = {};
        trialed{1,n}.time = {};
        trialed{1,n}.cfg = D_wavtransf{1,n}.cfg;
        trialed{1,n}.hdr = D_wavtransf{1,n}.hdr;
        trialed{1,n}.fsample = D_wavtransf{1,n}.fsample;
        %}
        
        % paste the data into the new table under the correct trial number
        % run through each trial
    %disp ('entering loop');
    for j=1:annot_sz(1) % trial
        %fprintf('j: %d', j);
        trial_temp = [];
        count = 1;

        trial_temp = D_wavtransf.trial{1,1}(:,timepoints(j,1):timepoints(j,2),:);
        time_temp = D_wavtransf.time{1,1}(1, timepoints(j,1):timepoints(j,2));

        %{
        for n = 1:numFreq
            %fprintf('frequency (n): %d', n);
            for k=timepoints(j,1):timepoints(j,2) % timepoints
                % take each datapoint between onset and offset from continous and put into .trial
                % for each electrode

                %cont_sz = size(continuous{1,n}.trial);
                cont_sz = size(D_wavtransf.trial{1,1}(:,:,n));
                for m=1:cont_sz(1);
                    %trial_temp(m,count) = continuous{1,n}.trial(m,k); % need to do for each electrode
                    trial_temp(m,count,n) = D_wavtransf.trial{1,1}(m,k,n); 
                end
                % take each timepoint between onset and offset and put into .time
                %time_temp(1,count) = continuous{1,n}.time(1,k);
                time_temp(1,count) = D_wavtransf.time{1,1}(1,k);
                count = count+1;
            end
        end
        %}
        %trialed{1,n}.trial{1,j} = trial_temp;
        %trialed{1,n}.time{1,j} = time_temp;

        trialed.trial{1,j} = trial_temp;
        trialed.time{1,j} = time_temp;
    end

    save([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_ft-raw_freqbands-trial.mat'],'trialed','-v7.3');
    %save([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_test_change-frequency-location'],'trialed','-v7.3');
    %end