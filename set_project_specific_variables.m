% this script should be called by functions within ieeg_ft_funcs_am....
% ... it will determine which data collection site (Pitt vs. MGH) the subject is from
% ... and then set paths and other variables appropriate to that site
%%%% implemented as script rather than function so that these vars will be put directly into workspace 

%% a lot of the paths that get set below can be replaced by calling setpaths_dbs_triplet and setpaths_dbs_seq

vardefault('op',struct);
field_default('op','sub','DBS3012'); % triplet subject
field_default('op','art_crit','E'); % default to high gamma 
proj_str = regexprep(op.sub, '[0-9]', ''); % string that will tell us what project / collection site this subject is from

switch op.art_crit
    case 'E'
        op.resp_signal = 'hg'; % high gamma
    case 'F'
        op.resp_signal = 'beta';
end

if op.denoised
    op.denoise_string = '_denoised';
else
    op.denoise_string = '_not_denoised';
end

switch proj_str
    case 'DBS' % Pitt
        setpaths_dbs_triplet()

        TASK = 'triplet'; 
        PATH_SUBJECT=[PATH_DATA filesep op.sub];
        PATH_DER_SUB = [PATH_SUBJECT filesep 'Preprocessed Data']; % subject-specific derivatives
        PATH_SYNC = [PATH_DER_SUB filesep 'Sync'];
        PATH_ANNOT = [PATH_SYNC '/annot']; 
        PATH_FIELDTRIP = [PATH_DER_SUB filesep 'FieldTrip'];
        PATH_FIGURES = [PATH_SYNC, filesep, 'figures']; 
        
        % filepath and filename for saving artifact info into subject-specific annot folder
        ARTIFACT_FILENAME_SUB = [PATH_ANNOT, filesep, op.sub, '_artifact_criteria_', op.art_crit, op.denoise_string, '.txt']; 
       
        % string (including filepath) at beginning of all fieldtrip filenames for this subj
        FT_FILE_PREFIX = [PATH_FIELDTRIP, filesep, op.sub, '_ft_'];

        FT_RAW_FILENAME = [FT_FILE_PREFIX 'raw_session.mat']; 
      
        artparam = readtable([PATH_ARTIFACT, filesep, 'artifact_', op.art_crit, '_params.txt']);
        session= bml_annot_read([PATH_ANNOT filesep op.sub '_session.txt']);
        electrodes = bml_annot_read([PATH_ANNOT filesep op.sub '_electrode.txt']);
            electrodes.name = electrodes.electrode; % match the table variable name used in dbs-seq
     
        %%% trial_epoch variable might only be used in analyzing vibration-denoised data
        % DBS3031 is missing DBS3031_trial_epoch.txt in annot folder; instead only has DBS3031_trial_epoch_criteria_D.txt
        if ~strcmp(op.sub, 'DBS3031') 
            trial_epoch = readtable([PATH_ANNOT filesep op.sub '_trial_epoch.txt']);
        elseif strcmp(op.sub, 'DBS3031') 
            trial_epoch = readtable([PATH_ANNOT filesep op.sub '_trial_epoch_criteria_D.txt']);
        end

        load_triplet_stim_beh_timing()

    case 'DM' % MGH
        setpaths_dbs_seq()
    
        SESSION = 'intraop';
        TASK = 'smsl'; 
        
        PATH_DER_SUB = [PATH_DER filesep 'sub-' op.sub];  
        PATH_PREPROC = [PATH_DER_SUB filesep 'preproc'];
        PATH_ANNOT = [PATH_DER_SUB filesep 'annot'];
        PATH_FIELDTRIP = [PATH_DER_SUB filesep 'fieldtrip'];
        PATH_AEC = [PATH_DER_SUB filesep 'aec']; 
        PATH_SCORING = [PATH_DER_SUB filesep 'analysis' filesep 'task-', TASK, '_scoring'];
        PATH_ANALYSIS = [PATH_DER_SUB filesep 'analysis'];
        PATH_TRIAL_AUDIO = [PATH_ANALYSIS filesep 'task-', TASK, '_trial-audio'];
        PATH_TRIAL_AUDIO_INTRAOP_GO = [PATH_TRIAL_AUDIO filesep 'ses-', SESSION, '_go-trials'];
        PATH_TRIAL_AUDIO_INTRAOP_STOP = [PATH_TRIAL_AUDIO filesep 'ses-', SESSION, '_stop-trials']; 
        
        PATH_SRC_SUB = [PATH_SRC filesep 'sub-' op.sub];  
        PATH_SRC_SESS = [PATH_SRC_SUB filesep 'ses-' SESSION]; 
        PATH_AUDIO = [PATH_SRC_SESS filesep 'audio']; 
        PATHS_TASK = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'task');

        % filepath and filename for saving artifact info into subject-specific annot folder
        ARTIFACT_FILENAME_SUB = [PATH_ANNOT filesep 'sub-' op.sub '_ses-' SESSION '_task-' TASK '_artifact-criteria-' op.art_crit, op.denoise_string, '.tsv'];

        % string (including filepath) at beginning of all fieldtrip filenames for this subj
        FT_FILE_PREFIX = [PATH_FIELDTRIP, filesep, 'sub-', op.sub, '_ses-' SESSION '_task-' TASK, '_ft-'];
        
        FT_RAW_FILENAME = [FT_FILE_PREFIX 'raw.mat']; 
        
        PATH_ART_PROTOCOL = ['Y:\DBS\groupanalyses\task-smsl\A09_artifact_criteria_', op.art_crit]; 
        PATH_FIGURES = [PATH_ART_PROTOCOL filesep 'figures']; 
        
        artparam = bml_annot_read_tsv([PATH_ARTIFACT, filesep, 'artifact_', op.art_crit, '_params.tsv']); 
        session= bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' op.sub '_sessions.tsv']);
        
        % merge electrode info into one table
        electrodes = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' op.sub '_electrodes.tsv']);
        channels = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' op.sub '_ses-' SESSION '_channels.tsv']); %%%% for connector info
            channels.name = strrep(channels.name,'_Ll','_Lm'); % change name to match naming convention in electrodes table
        [~, ch_ind] = intersect(channels.name, electrodes.name,'stable');
        electrodes = join(electrodes,channels(ch_ind,{'name','connector'}) ,'keys','name'); %%% add connector info

        % define trial epochs for referencing
        % for dbs-seq/smsl, we will use experimenter keypress for trial start/end times
        %%%%% this means no trial overlap, but generally a large time buffer before cue onset and after speech offset
        epoch = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' op.sub '_ses-' SESSION '_task-' TASK '_annot-trials.tsv']);
    otherwise
        error('Could not identify project from subject name')
end

% common variables
if exist(ARTIFACT_FILENAME_SUB, 'file')
    artifact = bml_annot_read(ARTIFACT_FILENAME_SUB);  
end



clear proj_str