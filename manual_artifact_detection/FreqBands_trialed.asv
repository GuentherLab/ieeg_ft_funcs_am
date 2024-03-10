SUBJECT = 'sub-DM1005';

if ~isempty(regexp(SUBJECT,'DBS30\d\d')) % if SUBJECT is triplet
    disp('triplet subject');

    PATH_DATASET = 'Z:\DBS';
    PATH_SUB = [PATH_DATASET filesep SUBJECT];
    PATH_PREPROC_DATA = [PATH_SUB filesep 'Preprocessed Data'];
    PATH_FT = [PATH_PREPROC_DATA filesep 'FieldTrip'];
    PATH_SYNC = [PATH_PREPROC_DATA filesep 'Sync'];
    PATH_ANNOT = [PATH_SYNC filesep 'annot'];
    
    filename = [SUBJECT '_ft_raw_session_freqbands'];
    
    % load freqband file
    tempname = load([PATH_FT filesep filename]);
    FT_file_fb = tempname.D_wavtransf;
    
    numFreq = length(FT_file_fb);
    
    % load annot file
    annot = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch'],'Delimiter','\t');
    
    % create new table with freqband file.trial and .time stored continuously
    continuous = cell(1,numFreq);
    for i = 1:numFreq
        continuous{1,i} = struct('time','trial');
        continuous{1,i}.time = horzcat(FT_file_fb{1,i}.time{1,1}, FT_file_fb{1,i}.time{1,2}, ...
            FT_file_fb{1,i}.time{1,3}, FT_file_fb{1,i}.time{1,4});
        continuous{1,i}.trial = horzcat(FT_file_fb{1,i}.trial{1,1}, FT_file_fb{1,i}.trial{1,2}, ...
            FT_file_fb{1,i}.trial{1,3}, FT_file_fb{1,i}.trial{1,4});
    end
    
    disp('calculating timepoints');
    % find all timepoints between start and end of each trial
    timepoints = []; % each row is a trial; first column is start location; second column is end location
        % both timepoints are from continuous
    annot_sz = size(annot);
    j=1;
    for i = 1:annot_sz(1)
        % for each timepoint (a) iterate through continuous and for each column in continuous (b) do a-b
        % when a-b becomes < 0, the previous timepoint is the selected column
        % optimization: after one timepoint is calculated, the calculation does not need to be performed on the previous set of times in continuous
    
        % loop to determine onset column
        onset=1;
        while onset<annot{i,2}
            onset = continuous{1,1}.time(1,j);
            j=j+1;
        end
        timepoints(i,1) = j;
        
        offset=1;
        while offset<annot{i,3}
            offset = continuous{1,1}.time(1,j);
            j=j+1;
        end
        timepoints(i,2) = j;
    end
    
    trialed = cell(1,numFreq);
    % loop through each frequency
    for i = 1:numFreq
        % copy other data stored in the table to the new variable (basically not .time and .trial)
        trialed{1,i} = struct('label','trial','time','cfg','hdr','fsample');
        trialed{1,i}.label = FT_file_fb{1,i}.label;
        trialed{1,i}.trial = {};
        trialed{1,i}.time = {};
        trialed{1,i}.cfg = FT_file_fb{1,i}.cfg;
        trialed{1,i}.hdr = FT_file_fb{1,i}.hdr;
        trialed{1,i}.fsample = FT_file_fb{1,i}.fsample;
        
        % paste the data into the new table under the correct trial number
        % run through each trial
        for j=1:annot_sz(1) % trial
            trial_temp = [];
            count = 1;
            for k=timepoints(j,1):timepoints(j,2) % timepoints
                % take each datapoint between onset and offset from continous and put into .trial
                % for each electrode
                cont_sz = size(continuous{1,i}.trial);
                for m=1:cont_sz(1);
                    trial_temp(m,count) = continuous{1,i}.trial(m,k); % need to do for each electrode
                end
                
                % take each timepoint between onset and offset and put into .time
                time_temp(1,count) = continuous{1,i}.trial(1,k);
                count = count+1;
            end
            trialed{1,i}.trial{1,j} = trial_temp;
            trialed{1,i}.time{1,j} = time_temp;
        end
    end
    
    save([PATH_FT filesep SUBJECT '_ft_raw_session_freqbands_trialed.mat'],'trialed','-v7.3');
elseif ~isempty(regexp(SUBJECT,'sub-DM10\d\d')) % else if SUBJECT is SMSL
    disp('SMSL subject');

    PATH_DATASET = 'Y:\DBS';
    PATH_DER = [PATH_DATASET filesep 'derivatives'];
    PATH_SUB = [PATH_DER filesep SUBJECT];
    PATH_FT = [PATH_SUB filesep 'fieldtrip'];
    PATH_ANNOT = [PATH_SUB filesep 'annot'];
    
    filename = [SUBJECT '_ses-intraop_task-smsl_ft-raw_freqbands'];
    
    % load freqband file
    tempname = load([PATH_FT filesep filename]);
    FT_file_fb = tempname.D_wavtransf;
    
    numFreq = length(FT_file_fb);
    
    % load annot file
    annot = readtable([PATH_ANNOT filesep SUBJECT '_ses-intraop_task-smsl_annot-trials.tsv'],'FileType','text','Delimiter','\t');
    
    % create new table with freqband file.trial and .time stored continuously
    continuous = cell(1,numFreq);
    for i = 1:numFreq
        continuous{1,i} = struct('time','trial');
        continuous{1,i}.time = FT_file_fb{1,i}.time{:,:};
        continuous{1,i}.trial = FT_file_fb{1,i}.trial{:,:};
    end

    disp('calculating timepoints');
    % find all timepoints between start and end of each trial
    timepoints = []; % each row is a trial; first column is start location; second column is end location
        % both timepoints are from continuous
    annot_sz = size(annot);
    j=1;
    for i = 1:annot_sz(1)
        % for each timepoint (a) iterate through continuous and for each column in continuous (b) do a-b
        % when a-b becomes < 0, the previous timepoint is the selected column
        % optimization: after one timepoint is calculated, the calculation does not need to be performed on the previous set of times in continuous
    
        % loop to determine onset column
        onset=1;
        while onset<annot{i,1}
            onset = continuous{1,1}.time(1,j);
            j=j+1;
        end
        timepoints(i,1) = j;
        
        offset=1;
        offsettimep = annot{i,1} + annot{i,2}; % annot{i,2} is the duration, not offset
        while offset<offsettimep
            offset = continuous{1,1}.time(1,j);
            j=j+1;
        end
        timepoints(i,2) = j;
    end
    
    trialed = cell(1,numFreq);
    % loop through each frequency
    for i = 1:numFreq
        % copy other data stored in the table to the new variable (basically not .time and .trial)
        trialed{1,i} = struct('label','trial','time','cfg','hdr','fsample');
        trialed{1,i}.label = FT_file_fb{1,i}.label;
        trialed{1,i}.trial = {};
        trialed{1,i}.time = {};
        trialed{1,i}.cfg = FT_file_fb{1,i}.cfg;
        trialed{1,i}.hdr = FT_file_fb{1,i}.hdr;
        trialed{1,i}.fsample = FT_file_fb{1,i}.fsample;
        
        % paste the data into the new table under the correct trial number
        % run through each trial
        for j=1:annot_sz(1) % trial
            trial_temp = [];
            count = 1;
            for k=timepoints(j,1):timepoints(j,2) % timepoints
                % take each datapoint between onset and offset from continous and put into .trial
                % for each electrode
                cont_sz = size(continuous{1,i}.trial);
                for m=1:cont_sz(1);
                    trial_temp(m,count) = continuous{1,i}.trial(m,k); % need to do for each electrode
                end
                
                % take each timepoint between onset and offset and put into .time
                time_temp(1,count) = continuous{1,i}.trial(1,k);
                count = count+1;
            end
            trialed{1,i}.trial{1,j} = trial_temp;
            trialed{1,i}.time{1,j} = time_temp;
        end
    end
    
    save([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_ft-raw_freqbands_trialed.mat'],'trialed','-v7.3');
else
    disp('invalid subject');
end

