%% SMSL
SUBJECT = 'sub-DM1005';

PATH_DATASET = 'Y:\DBS';
PATH_DER = [PATH_DATASET filesep 'derivatives'];
PATH_SUB = [PATH_DER filesep SUBJECT];
PATH_FT = [PATH_SUB filesep 'fieldtrip'];
PATH_ANNOT = [PATH_SUB filesep 'annot'];

filename = [SUBJECT '_ses-intraop_task-smsl_ft-raw'];

tempname = load([PATH_FT filesep filename]);
FT_file = tempname.D;

annot = readtable([PATH_ANNOT filesep SUBJECT '_ses-intraop_task-smsl_annot-produced-syllables.tsv'],'FileType','text','Delimiter','\t');

annot_sz = size(annot);
epoch_table = table('Size',[annot_sz(1) 2],'VariableTypes',["double","double"]);
epoch_table.Properties.VariableNames = ["starts","ends"];
epoch_table.starts(:) = annot.audio_onset(:) - 0.5;
for i=1:length(annot.sp_off)
    if isnan(annot.sp_off(i)) && ~isnan(annot.keypress_time(i))
        epoch_table.ends(i) = annot.keypress_time(i); % if the trial is unusable and there is a keypress time it will instead get data ending at the keypress time
    elseif isnan(annot.sp_off(i)) && isnan(annot.keypress_time(i))
        epoch_table.ends(i) = (annot.audio_onset(i) - 0.5) + 5;
    else
        epoch_table.ends(i) = annot.sp_off(i) + 1.5;
    end
end

cfg = [];
cfg.epoch = epoch_table;
D_trialed = bml_redefinetrial(cfg, FT_file);

save([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_ft-raw_trial'],'D_trialed','-v7.3');

%% triplet
SUBJECT = 'DBS3001';
disp(SUBJECT);

clear trial_startend
clear D_trialed

PATH_DATASET = 'Z:\DBS';
PATH_SUB = [PATH_DATASET filesep SUBJECT];
PATH_PREPROC_DATA = [PATH_SUB filesep 'Preprocessed Data'];
PATH_FT = [PATH_PREPROC_DATA filesep 'FieldTrip'];
PATH_SYNC = [PATH_PREPROC_DATA filesep 'Sync'];
PATH_ANNOT = [PATH_SYNC filesep 'annot'];

filename = [SUBJECT '_ft_raw_session.mat'];

tempname = load([PATH_FT filesep filename]);
FT_file = tempname.D;
FT_file.time = FT_file.time.';
FT_file.trial = FT_file.trial.';

annot = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch'],'Delimiter','\t');

% make the data continuous (originally stored in 4 groups)
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

annot_syllable = readtable([PATH_ANNOT filesep SUBJECT '_produced_syllable'],'Delimiter','\t');

% create a list of only the end of the last syllable
annot_end = annot_syllable.ends(3:3:end);

% where should I get the trial window from for triplet subjects?
annot_stim = bml_annot_read_tsv([PATH_ANNOT filesep SUBJECT '_stimulus_triplet']);

disp('determining where the trials start and end');
% determine where the trials start and end
%annot_sz = size(annot);
%epoch_table = table('Size',[annot_sz(1) 2],'VariableTypes',["double","double"]);
%epoch_table.Properties.VariableNames = ["starts","ends"];
%epoch_table.starts(:) = annot_stim.starts(:) - 0.5;

%{
for i=1:length(annot_end)
    rownum_thirdsyl = i*3; % row number for the third syllable in _produced_syllable

    if ~isnan(annot_end(i)) % third syllable offset
        % if the speech off of the third syllable exists, use that
        epoch_table.ends(i) = annot_end(i) + 1.5;
    elseif ~isnan(annot_syllable.starts(rownum_thirdsyl)) % third syllable onset
        epoch_table.ends(i) = annot_syllable.starts(rownum_thirdsyl) + 1.5;
    elseif ~isnan(annot_syllable.ends(rownum_thirdsyl-1)) % second syllable offset
        epoch_table.ends(i) = annot_syllable.ends(rownum_thirdsyl-1) + 1.5;
    elseif ~isnan(annot_syllable.starts(rownum_thirdsyl-1)) % second syllable onset
        epoch_table.ends(i) = annot_syllable.starts(rownum_thirdsyl-1) + 1.5;
    elseif ~isnan(annot_syllable.ends(rownum_thirdsyl-2)) % first syllable offset
        epoch_table.ends(i) = annot_syllable.ends(rownum_thirdsyl-2) + 1.5;
    elseif ~ isnan(annot_syllable.starts(rownum_thirdsyl-2)) % first syllable onset
        epoch_table.ends(i) = annot_syllable.starts(rownum_thirdsyl-2) + 1.5;
    else
        fprintf('no triplet speech onset and offset data, used stim offset. trial number: %d \n', i);
        epoch_table.ends(i) = annot_stim.ends(i) + 4;
    end
end
%}

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

D_trialed.trial = cell(length(FT_file.trial),length(trial_startend));
D_trialed.time = cell(length(FT_file.trial),length(trial_startend));

disp('cutting up into trials');
for i=1:length(FT_file.trial) % session
    temp_startend = trial_startend(:,:,i);
    delRows = find(trial_startend(:,1,i) == 0 | trial_startend(:,2,i) == 0);
    temp_startend(delRows,:) = [];

    epoch_table = table('Size',[length(temp_startend), 2],'VariableTypes',["double","double"]);
    epoch_table.Properties.VariableNames = ["starts","ends"];

    epoch_table.starts = temp_startend(:,1);
    epoch_table.ends = temp_startend(:,2);

    temp_FT.label = FT_file.label;
    temp_FT.hdr = FT_file.hdr;
    temp_FT.cfg = FT_file.cfg;
    temp_FT.trial = FT_file.trial(i,1);
    temp_FT.time = FT_file.time(i,1);

    cfg = [];
    cfg.epoch = epoch_table;
    temp = bml_redefinetrial(cfg,temp_FT); 

    D_trialed.trial{i,:} = temp.trial(1,:);
    D_trialed.time{i,:} = temp.time(1,:};
end

D_trialed.label = temp.label;
D_trialed.hdr = temp.hdr;
D_trialed.fsample = temp.fsample;

%save([PATH_FT filesep SUBJECT '_ft_raw_session_trial'],'D_trialed','-v7.3');
disp('finished');
