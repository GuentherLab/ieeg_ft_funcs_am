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

PATH_DATASET = 'Z:\DBS';
PATH_SUB = [PATH_DATASET filesep SUBJECT];
PATH_PREPROC_DATA = [PATH_SUB filesep 'Preprocessed Data'];
PATH_FT = [PATH_PREPROC_DATA filesep 'FieldTrip'];
PATH_SYNC = [PATH_PREPROC_DATA filesep 'Sync'];
PATH_ANNOT = [PATH_SYNC filesep 'annot'];

filename = [SUBJECT '_ft_raw_session.mat'];

tempname = load([PATH_FT filesep filename]);
FT_file = tempname.D;

annot = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch'],'Delimiter','\t');

% make the data continuous (originally stored in 4 groups)
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


