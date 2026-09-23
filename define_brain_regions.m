%%%% assign broad region definitions based on small area labels
% 
% op.atlas must not be empty; for MGH BML projects this will generally be hcp_distal
%   have not added definitions here yet for triplet/pitt

% if optional labels_in table of electrodes is provided, then region labels will be added in table variable 'region'
%%% this labels_in table will also be used to figure out which atlas we are using

function [labels_out, op_out] = define_brain_regions(labels_in, op)

vardefault('labels_in', {}); 
vardefault('op', struct); 
field_default('op', 'include_bottom_all_row', 0); % add bottom row labeled 'all'

% process labels_in variable type
argtype = class(labels_in);
if isempty(labels_in) || (iscell(labels_in) && numel(labels_in) == 0)
    argtype = 'empty';
end

labels_out = labels_in;

% check if we're using the HCP atlas 
if strcmp(argtype, 'table')
    if any(contains(labels_in.Properties.VariableNames, 'HCPMMP1_label_1')) % this variable should usually be in labels_in table for MGH data
        op.atlas = 'hcp_distal'; 
    end
end
field_default('op', 'atlas', 'hcp_distal') % use atlas default if not yet defined

% Define region-area mappings based on atlas
switch op.atlas 
    case 'hcp_distal'
        % 1 = 'Area 1' (Fischl et al 2008, Geyer et al 1999, Geyer et al 2000) ... posterior postcentral gyrus
        % 2 = 'Area 2' ... postcentral sulcus
        % 3a = 'Area 3a'.... central sulcus
        % 3b = 'Primary Sensory Cortex'.... postcentral gyrus
        % 4 = 'Primary Motor Cortex'.... anterior central sulcus
        % 6v = 'Ventral Area 6' (Fischl et al 2008, Amunts et al 2010, Geyer 2004).... precentral gyrus, precentral sulcus
        %
        % 6d = 'Dorsal Area 6' (Fischl et al 2008, Geyer 2004, Geyer et al 2000).... dorsal precentral gyrus (hand knob?)
        % 
        % 6r = 'Rostral Area 6' (Amunts et al 2010)... ventral premotor, precentral sulcus
        %
        % FEF = Frontal Eye Fields... in first 5 dbsseq subs, this is close to precentral gyrus, but may be more frontal in future subs
        % PEF = Premotor Eye Fields
        %
        % OP4 = 'Area OP4/PV' .... ventral precentral/postecentral gyrus, operculum
        %
        % 55b = 'Area 55b' (Hopf 1956)... mid precentral gyrus, precentral sulcus, posterior MFG... premotor cortex
        % 
        % 43 = 'Area 43' (Brodmann 1909, Brodmann 2007, Nieuwenhuys et al 2014)... operculum and ventral precentral gyrus
        %
        % i6-8 = 'Inferior 6-8 Transitional Area'(von Economo and Koskinas 1925, Triarhou 2007)... dorsal premotor
        %
        % 8Av = 'Area 8av' (Petredes and Pandya 1999) .... middle frontal gyrus
        % 8C = 'Area 8C' (Petredes and Pandya 1999) ... ventral middle frontal gyrus
        %
        % A4 = 'Auditory 4 Complex' (Morosan et al 2005).... dorsal STG
        % A5 = 'Auditory 5 Complex' .... ventral STG
        %
        % PF = 'Area PF Complex'.... supramarginal gyrus
        % PFop = 'Area PF opercular'... ant supramarginal gyrus, operculum, ventral postcentral sulcus
        %
        % PSL = 'PeriSylvian Language Area'.... angular gyrus
        %
        % TE1a = 'Area TE1 anterior' (von Economo and Koskinas 1925, Triarhou 2007)... ant middle temporal gyrus
        %
        % STV = 'Superior Temporal Visual Area' .... post STG

        region_areas = {   'SMC',  {'1','2','3a','3b','4','6v','6d','43','55b','PEF','FEF','OP4','i6-8'};...
                        'STG', {'A4','A5','STGa','STV','TPOJ1'};...
                        'MFG',  {'8Av','8C','p9-46v'};...
                        'IFG/IFS',  {'44','45','IFSp'};...
                        'SMG/PF', {'PF','PFop'};...
                        'MTG', {'TE1a','TE1m','TE1p'};...
                        'STN', {'STN_associative_L','STN_motor_L','STN_motor_R'};...
                        'Thal', {'087_Thalamus_ventro_oralis_anterior_Voa_L','088_Thalamus_ventro_oralis_posterior_Vop_L','088_Thalamus_ventro_oralis_posterior_Vop_R',...
                                   '090_Thalamus_zentrolateralis_oralis_Zo_L','091_Thalamus_ventro_intermedius_internus_Vimi_R',...
                                   '094_Thalamus_ventro_intermedius_externus_Vime_L','094_Thalamus_ventro_intermedius_externus_Vime_R'};...
                        'GP', {'GPe_L','GPe_R','GPi_postparietal_R','GPi_premotor_R','GPi_sensorimotor_L','GPi_sensorimotor_R'};...
                        };

        op_out.atlas_var_names = {'HCPMMP1_label_1';'DISTAL_label_1'}; 

    otherwise  % no other atlases implemented yet
        error('unrecognized atlas') 

end

% Initialize regiondef table with electode lists for each region
op_out.regiondef = table(region_areas(:,1), region_areas(:,2), cell(size(region_areas,1),1), ...
                            'VariableNames', {'region', 'areas', 'electrode_list'});
op_out.nregions = height(op_out.regiondef); 

% Process labels based on input type
switch argtype 
    case 'empty'
        % No input provided, just return the region definitions
        labels_out = {};
        
    case 'cell' 
        % Input is cell array of area labels
        % Output is cell array of region labels (same size as input, with indices tracked)
        labels_out = cell(size(labels_in));
        
        for ielc = 1:length(labels_in)
            area_label = labels_in{ielc};
            % Find which region this area belongs to
            for iregion = 1:op_out.nregions
                if any(strcmp(area_label, op_out.regiondef.areas{iregion}))
                    labels_out{ielc} = op_out.regiondef.region{iregion};
                    % Track this electrode index in the region
                    op_out.regiondef.electrode_list{iregion} = [op_out.regiondef.electrode_list{iregion}; ielc];
                    break;
                end
            end
            if isempty(labels_out{ielc})
                labels_out{ielc} = 'unknown'; % area not found in any region
            end
        end
        
    case 'table' 
        % Expected to have at least one of the op_out.atlas_var_names as a table variable
        labels_out = labels_in; 
        labels_out.region = cell(height(labels_out), 1);
        
        % Determine which column contains electrode names (if available)
        if any(contains(labels_out.Properties.VariableNames, 'chan'))
            elc_name_var = 'chan';
        elseif any(contains(labels_out.Properties.VariableNames, 'name'))
            elc_name_var = 'name';
        else
            elc_name_var = []; % no electrode name column found
        end
        
        % Move region column to the front (before atlas variables)
        labels_out = movevars(labels_out, 'region', 'Before', op_out.atlas_var_names{1}); 
        
        % Assign regions based on area labels
        for iregion = 1:op_out.nregions
            thisregion = op_out.regiondef.region{iregion};
            
            for iatlas = 1:length(op_out.atlas_var_names) 
                atlas_var = op_out.atlas_var_names{iatlas};
                
                % Check if this atlas variable exists in the table
                if any(strcmp(labels_out.Properties.VariableNames, atlas_var))
                    elcs_in_this_region = ismember(labels_out{:, atlas_var}, op_out.regiondef.areas{iregion}); 
                    labels_out.region(elcs_in_this_region) = {thisregion};
                    
                    % Track electrode names/indices in this region
                    if ~isempty(elc_name_var)
                        op_out.regiondef.electrode_list{iregion} = ...
                            [op_out.regiondef.electrode_list{iregion}; labels_out{elcs_in_this_region, elc_name_var}];
                    else
                        op_out.regiondef.electrode_list{iregion} = ...
                            [op_out.regiondef.electrode_list{iregion}; find(elcs_in_this_region)];
                    end
                end
            end
        end
        
    otherwise 
        error('unrecognized electrode list type')
end

% Optionally add bottom row for convenience, for calling functions to use when compiling analyses across regions
if op.include_bottom_all_row
    regiondef_varnames = op_out.regiondef.Properties.VariableNames; 
    all_elcs = [];
    for iregion = 1:op_out.nregions
        all_elcs = [all_elcs; op_out.regiondef.electrode_list{iregion}];
    end
    op_out.regiondef = [op_out.regiondef; table({'all'}, {{}}, {all_elcs}, 'VariableNames', regiondef_varnames, 'RowNames', {'all'})]; 
end

% Copy op to op_out (preserve any fields from input op)
op_out = mergestruct(op_out, op);

end

function out = mergestruct(varargin)
% Helper function to merge structures, with later arguments taking precedence
out = varargin{1};
for i = 2:nargin
    fnames = fieldnames(varargin{i});
    for j = 1:length(fnames)
        out.(fnames{j}) = varargin{i}.(fnames{j});
    end
end
end