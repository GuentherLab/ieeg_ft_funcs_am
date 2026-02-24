%%%% Andrew Meier matlab startup script for analysis of MGH BML data on Turbo server

% map network drives
!net use Y: \\132.183.240.49\Nexus4
!net use Z: \\132.183.240.28\Nexus

%%% add paths
full_paths_to_add = {'C:\Program Files\Brain-Modulation-Lab\bml';.... % old location was 'Y:\Documents\Code\bml';
                'Y:\Documents\Code\dbs_seq_analysis';...
                'C:\Program Files\LeadDBS_Classic';... % use this rather than updated LeadDBS for MGH preprocessing
%                 'C:\Users\amsmeier\dbs_triplet'...
                };
single_paths_to_add = {'C:\Program Files\NPMK';...
                        'Y:\Users\lbullock\MATLAB_external_libs_Turbo20230907\fieldtrip';... % latane's copy of fieltrip; specific version... do not use most up-to-date version
                        'Y:\DBS\derivatives\sub-DM-example\preproc\utils';... % used in MGH preproc step B06
                        'C:\Program Files\SPM\spm12';...
                        }; 
full_paths_to_remove = {'Y:\Documents\Code\dbs_seq_analysis\archive';...
                        'C:\Users\amsmeier\dbs_triplet\archive\preprocessing';... 
                        'C:\Users\amsmeier\dbs_triplet\archive';...
                        'C:\Program Files\LeadDBS'}; 

for ipath = 1:length(full_paths_to_add)
    addpath(genpath(full_paths_to_add{ipath}))
end
for ipath = 1:length(single_paths_to_add)
    addpath(single_paths_to_add{ipath})
end
for ipath = 1:length(full_paths_to_remove)
    rmpath(genpath(full_paths_to_remove{ipath}))
end


ft_defaults() % function from top level of fieldtrip
bml_defaults()

%%%%% go to main dir
cd('Y:\DBS\derivatives')
% cd('Z:\DBS')

clear