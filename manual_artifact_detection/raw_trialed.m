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
epoch_table = zeros([annot_sz(1),2]);
epoch_table(:,1) = annot.audio_onset(:) - 0.5;
for i=1:length(annot.sp_off)
    if isnan(annot.sp_off(i)) && ~isnan(annot.keypress_time(i))
        epoch_table(i,2) = annot.keypress_time(i); % if the trial is unusable and there is a keypress time it will instead get data ending at the keypress time
    elseif isnan(annot.sp_off(i)) && isnan(annot.keypress_time(i))
        epoch_table(i,2) = (annot.audio_onset(i) - 0.5) + 5;
    else
        epoch_table(i,2) = annot.sp_off(i) + 1.5;
    end
end

%{
cfg = [];
cfg.epoch = epoch_table;
D_trialed = bml_redefinetrial(cfg, FT_file);
%}

% calculating timepoints (column number)
timepoints = [];
j=1;
for k = 1:annot_sz(1)
    % for each timepoint (a) iterate through continuous and for each column in continuous (b) do a-b
    % when a-b becomes < 0, the previous timepoint is the selected column
    % optimization: after one timepoint is calculated, the calculation does not need to be performed on the previous set of times in continuous

    % loop to determine onset column
    onset=1;
    while onset<epoch_table(k,1)
        %onset = continuous{1,1}.time(1,j);
        onset = FT_file.time{1,1}(1,j);
        j=j+1;
    end
    timepoints(k,1) = j;
    
    stop=1;
    %offsettimep = annot{k,1} + annot{k,2}; % annot{i,2} is the duration, not offset
    while stop<epoch_table(k,2)
        %stop = continuous{1,1}.time(1,j);
        stop = FT_file.time{1,1}(1,j);
        j=j+1;
    end
    timepoints(k,2) = j;

    j = 1; % the trials may overlap, so it is necessary to reset j
end

% breaking up into trials
trialed = struct('label','trial','time','cfg','hdr','fsample');
trialed.label = FT_file.label;
trialed.trial = {};
trialed.time = {};
trialed.cfg = FT_file.cfg;
trialed.hdr = FT_file.hdr;
trialed.fsample = NaN;

for j=1:annot_sz(1) % trial
    %fprintf('j: %d', j);
    trial_temp = [];
    count = 1;
   
    trialed.trial{1,j} = FT_file.trial{1,1}(:,timepoints(j,1):timepoints(j,2));
    trialed.time{1,j} = FT_file.time{1,1}(1, timepoints(j,1):timepoints(j,2));
end

save([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_ft-raw_trial'],'trialed','-v7.3');
