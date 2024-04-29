triplet_subjects = [];
SMSL_subjects = ["sub-DM1005", "sub-DM1007", "sub-DM1008", "sub-DM1024", "sub-DM1025", "sub-DM1037"];

%% triplet loop
% need to change how the frequencies are stored (3d matrix in fieldtrip file .trial)
%for i=1:length(triplet_subjects)
    %SUBJECT = triplet_subjects(i); 
    SUBJECT = 'DBS3001';
    disp(SUBJECT);

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

    % size to trial matrix to determine size of transformed trial structure
    %{
    sz = size(continuous_FT.trial{1,1});

    envelope_wavpow_outfreq = 100;
    raw2wavpow_confactor = FT_file.fsample/envelope_wavpow_outfreq;
    new_sz = round(sz(2)/raw2wavpow_confactor);
    %}

    % calculate frequencies
    freqs = 10.^(linspace(0,2.5,20)); % number of frequencies: 20
    %D_wavtransf.trial{1,1} = zeros([sz(1),new_sz,length(freqs)]);
    for ifreq = 1:length(freqs)
        thisfreq = freqs(ifreq)
        cfg=[];
        cfg.out_freq = 100;
        cfg.wav_freq = thisfreq;
        cfg.wav_width = 7;

        %D_wavtransf{ifreq} = bml_envelope_wavpow(cfg,FT_file);
        temp = bml_envelope_wavpow(cfg,continuous_FT);

        D_wavtransf.trial{1,1}(:,:,ifreq) = temp.trial{1,1};
        D_wavtransf.time = temp.time;
        D_wavtransf.cfg = temp.cfg;
        D_wavtransf.cfg = rmfield(D_wavtransf.cfg,'wav_freq');
    end

    D_wavtransf.label = FT_file.label;
    D_wavtransf.hdr = FT_file.hdr;

    % break up into trials
    numFreq = length(D_wavtransf);
    annot = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch'],'Delimiter','\t');

    % save continuous (untrialed) data
    save([PATH_FT filesep SUBJECT '_ft_raw_session_freqbands.mat'],'trialed','-v7.3');

    % create new table with freqband file.trial and .time stored continuously
    %{
    continuous = cell(1,numFreq);
    for j = 1:numFreq
        continuous{1,j} = struct('time','trial');
        continuous{1,j}.time = horzcat(D_wavtransf{1,j}.time{1,1}, D_wavtransf{1,j}.time{1,2}, ...
            D_wavtransf{1,j}.time{1,3}, D_wavtransf{1,j}.time{1,4});
        continuous{1,j}.trial = horzcat(D_wavtransf{1,j}.trial{1,1}, D_wavtransf{1,j}.trial{1,2}, ...
            D_wavtransf{1,j}.trial{1,3}, D_wavtransf{1,j}.trial{1,4});
    end
    %}
    
    disp('calculating timepoints');
    % find all timepoints between start and end of each trial
    timepoints = []; % each row is a trial; first column is start location; second column is end location
        % both timepoints are from continuous
    annot_sz = size(annot);
    j=1;
    for k = 1:annot_sz(1)
        % for each timepoint (a) iterate through continuous and for each column in continuous (b) do a-b
        % when a-b becomes < 0, the previous timepoint is the selected column
        % optimization: after one timepoint is calculated, the calculation does not need to be performed on the previous set of times in continuous
    
        % loop to determine onset column
        onset=1;
        while onset<annot{k,2}
            onset = D_wavtransf.time{1,1}(1,j);
            j=j+1;
        end
        timepoints(k,1) = j;
        
        stop=1;
        while stop<annot{k,3}
            stop = D_wavtransf.time{1,1}(1,j);
            j=j+1;
        end
        timepoints(k,2) = j;
    end
    
    %trialed = cell(1,numFreq);
    % loop through each frequency
    %for n = 1:numFreq
        % copy other data stored in the table to the new variable (basically not .time and .trial)
        %trialed = struct('label','trial','time','cfg','hdr','fsample');
        trialed.label = D_wavtransf.label;
        trialed.trial = {};
        trialed.time = {};
        trialed.cfg = D_wavtransf.cfg;
        trialed.hdr = D_wavtransf.hdr;
        %trialed.fsample = D_wavtransf.fsample;
        
        % paste the data into the new table under the correct trial number
        % run through each trial
        for j=1:annot_sz(1) % trial
            trial_temp = [];
            count = 1;
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
            trialed.trial{1,j} = D_wavtransf.trial{1,1}(:,timepoints(j,1):timepoints(j,2),:);
            trialed.time{1,j} = D_wavtransf.time{1,1}(:,timepoints(j,1):timepoints(j,2));
        %end
    end

    save([PATH_FT filesep SUBJECT '_ft_raw_session_freqbands-trial.mat'],'trialed','-v7.3');
    %save([PATH_FT filesep SUBJECT 'ft_raw_session_test_change-frequency-location.mat'],'trialed','-v7.3');
%end

%% SMSL loop
%for n=1:length(SMSL_subjects)
    %SUBJECT = char(SMSL_subjects(i)); 
    SUBJECT = 'sub-DM1005';
    disp(SUBJECT);

    PATH_DATASET = 'Y:\DBS';
    PATH_DER = [PATH_DATASET filesep 'derivatives'];
    PATH_SUB = [PATH_DER filesep SUBJECT];
    PATH_FT = [PATH_SUB filesep 'fieldtrip'];
    PATH_ANNOT = [PATH_SUB filesep 'annot'];

    filename = [SUBJECT '_ses-intraop_task-smsl_ft-raw'];

    tempname = load([PATH_FT filesep filename]);
    FT_file = tempname.D;

    % computes trial start and end time
    %annot = readtable([PATH_ANNOT filesep SUBJECT '_ses-intraop_task-smsl_annot-trials.tsv'],'FileType','text','Delimiter','\t');
    annot = readtable([PATH_ANNOT filesep SUBJECT '_ses-intraop_task-smsl_annot-produced-syllables.tsv'],'FileType','text','Delimiter','\t');
    trial_startend(:,1) = annot.audio_onset(:) - 0.5;
    %trial_startend(:,2) = annot.sp_off(:) + 1.5;
    for i=1:length(annot.sp_off)
        if isnan(annot.sp_off(i)) && ~isnan(annot.keypress_time(i))
            trial_startend(i,2) = annot.keypress_time(i); % if the trial is unusable and there is a keypress time it will instead get data ending at the keypress time
        elseif isnan(annot.sp_off(i)) && isnan(annot.keypress_time(i))
            trial_startend(i,2) = (annot.audio_onset(i) - 0.5) + 5;
        else
            trial_startend(i,2) = annot.sp_off(i) + 1.5;
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

    freqs = 10.^(linspace(0,2.5,20)); % number of frequencies: 20
    D_wavtransf.trial{1,1} = zeros([sz(1),new_sz,length(freqs)]);
    for ifreq = 1:length(freqs)
        fprintf('ifreq: %d \n', ifreq);
        thisfreq = freqs(ifreq);
        cfg=[];
        cfg.out_freq = envelope_wavpow_outfreq;
        cfg.wav_freq = thisfreq;
        cfg.wav_width = 7;

        temp = bml_envelope_wavpow(cfg,FT_file);
        %D_wavtransf{ifreq} = bml_envelope_wavpow(cfg,FT_file); 

        D_wavtransf.trial{1,1}(:,:,ifreq) = temp.trial{1,1};
        D_wavtransf.time = temp.time;
        %{
        for iD = 1:sz(1)
            for jD = 1:sz(2)
                D_wavtransf.trial{1,1}(iD,jD,ifreq) = temp.trial{1,1}(iD,jD);
            end
        end
        %}
        D_wavtransf.cfg = temp.cfg;
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

