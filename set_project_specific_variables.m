% this script should be called by functions within ieeg_ft_funcs_am....
% ... it will determine which data collection site (Pitt vs. MGH) the subject is from
% ... and then set paths and other variables appropriate to that site
%%%% implemented as script rather than function so that these vars will be put directly into workspace 

%% a lot of the paths that get set below can be replaced by calling setpaths_dbs_triplet and setpaths_dbs_seq


vardefault('op',struct);
field_default('op','subj','DBS3012'); % triplet subject
field_default('op','art_crit','E'); % default to high gamma 
SUBJECT=op.sub;
proj_str = regexprep(SUBJECT, '[0-9]', ''); % string that will tell us what project / collection site this subject is from

switch proj_str
    case 'DBS' % Pitt
        setpaths_dbs_triplet()

        TASK = 'triplet'; 

        DATE=datestr(now,'yyyymmdd');
        PATH_DATA='Z:\DBS';
        PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
        PATH_DER_SUB = [PATH_SUBJECT filesep 'Preprocessed Data']; % subject-specific derivatives
        PATH_SYNC = [PATH_DER_SUB filesep 'Sync'];
        PATH_ART_PROTOCOL = ['Z:\DBS\Batch\P08_artifact_criteria_', op.art_crit];
        PATH_ANNOT = [PATH_SYNC '/annot']; 
        PATH_FIELDTRIP = [PATH_DER_SUB filesep 'FieldTrip'];
        
        % filepath and filename for saving artifact info into subject-specific annot folder
        ARTIFACT_FILENAME_SUB = [PATH_ANNOT filesep SUBJECT '_artifact_criteria_', op.art_crit, '_not-denoised.txt']; 
       
        % string (including filepath) at beginning of all fieldtrip filenames for this subj
        FT_FILE_PREFIX = [PATH_FIELDTRIP, filesep, SUBJECT, '_ft_'];

        FT_RAW_FILENAME = [FT_FILE_PREFIX 'raw_session.mat']; 
        
        artparam = readtable([PATH_ART_PROTOCOL, filesep, 'artifact_', op.art_crit , '_params.txt'],'FileType','text');
        session= bml_annot_read([PATH_ANNOT filesep SUBJECT '_session.txt']);
        electrodes = bml_annot_read([PATH_ANNOT filesep SUBJECT '_electrode.txt']);
            electrodes.name = electrodes.electrode; % match the table variable name used in dbs-seq
        
        %%% trial_epoch variable might only be used in analyzing vibration-denoised data
        % DBS3031 is missing DBS3031_trial_epoch.txt in annot folder; instead only has DBS3031_trial_epoch_criteria_D.txt
        if ~strcmp(SUBJECT, 'DBS3031') 
            trial_epoch = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch.txt']);
        elseif strcmp(SUBJECT, 'DBS3031') 
            trial_epoch = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch_criteria_D.txt']);
        end

        % define trial epochs for referencing
        % % % Use generous buffers around the cue and speech production, even if this means that consecutive trials overlap. 
        cue_presentation = bml_annot_read([PATH_ANNOT filesep SUBJECT '_cue_presentation.txt']);
        epoch = cue_presentation(:,{'stim1_onset','ends','session_id','trial_id'});
        epoch.starts = epoch.stim1_onset - 1.5;
        epoch.ends = epoch.ends + 2;
        epoch = bml_annot_table(epoch);

    case 'DM' % MGH
        setpaths_dbs_seq()
    
        SESSION = 'intraop';
        TASK = 'smsl'; 
        
        PATH_DER_SUB = [PATH_DER filesep 'sub-' SUBJECT];  
        PATH_PREPROC = [PATH_DER_SUB filesep 'preproc'];
        PATH_ANNOT = [PATH_DER_SUB filesep 'annot'];
        PATH_FIELDTRIP = [PATH_DER_SUB filesep 'fieldtrip'];
        PATH_AEC = [PATH_DER_SUB filesep 'aec']; 
        PATH_SCORING = [PATH_DER_SUB filesep 'analysis' filesep 'task-', TASK, '_scoring'];
        PATH_ANALYSIS = [PATH_DER_SUB filesep 'analysis'];
        PATH_TRIAL_AUDIO = [PATH_ANALYSIS filesep 'task-', TASK, '_trial-audio'];
        PATH_TRIAL_AUDIO_INTRAOP_GO = [PATH_TRIAL_AUDIO filesep 'ses-', SESSION, '_go-trials'];
        PATH_TRIAL_AUDIO_INTRAOP_STOP = [PATH_TRIAL_AUDIO filesep 'ses-', SESSION, '_stop-trials']; 
        
        PATH_SRC_SUB = [PATH_SRC filesep 'sub-' SUBJECT];  
        PATH_SRC_SESS = [PATH_SRC_SUB filesep 'ses-' SESSION]; 
        PATH_AUDIO = [PATH_SRC_SESS filesep 'audio']; 
        PATHS_TASK = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'task');

        % filepath and filename for saving artifact info into subject-specific annot folder
        ARTIFACT_FILENAME_SUB = [PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_task-' TASK '_artifact-criteria-' ARTIFACT_CRIT '_not-denoised.tsv'];

        % string (including filepath) at beginning of all fieldtrip filenames for this subj
        FT_FILE_PREFIX = [PATH_FIELDTRIP, filesep, 'sub-', SUBJECT, '_ses-' SESSION '_task-' TASK, '_ft-'];
        
        FT_RAW_FILENAME = [FT_FILE_PREFIX 'raw.mat']; 
        
        PATH_ART_PROTOCOL = ['Y:\DBS\groupanalyses\task-smsl\A09_artifact_criteria_', ARTIFACT_CRIT];
        PATH_FIGURES = [PATH_ART_PROTOCOL filesep 'figures']; 

        artparam = readtable([PATH_ART_PROTOCOL, filesep, 'artifact_', ARTIFACT_CRIT , '_params.tsv'],'FileType','text');
        session= bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_sessions.tsv']);

        % merge electrode info into one table
        electrodes = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_electrodes.tsv']);
        channels = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_channels.tsv']); %%%% for connector info
            channels.name = strrep(channels.name,'_Ll','_Lm'); % change name to match naming convention in electrodes table
        [~, ch_ind] = intersect(channels.name, electrodes.name,'stable');
        electrodes = join(electrodes,channels(ch_ind,{'name','connector'}) ,'keys','name'); %%% add connector info

        % define trial epochs for referencing
        % for dbs-seq/smsl, we will use experimenter keypress for trial start/end times
        %%%%% this means no trial overlap, but generally a large time buffer before cue onset and after speech offset
        epoch = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_task-' TASK '_annot-trials.tsv']);
    otherwise
        error('Could not identify project from subject name')
end

% common variables





clear proj_str