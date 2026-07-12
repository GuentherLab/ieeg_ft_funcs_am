%%%% assign broad region definitions based on small area labels
% 
% op.atlas must not be empty; for MGH BML projects this will generally be hcp_distal
%   have not added definitions here yet for triplet/pitt

% if optional 'resp' table of electrodes is provided, then region labels will be added in table variable 'region'
%%% this resp table will also be used to figure out which atlas we are using

function [op_out, resp_out] = define_brain_regions(op,resp)

op_out = op; 

if exist('resp','var')
    if any(contains(resp.Properties.VariableNames,'HCPMMP1_label_1')) % this variable should usually be in resp table for MGH data
        op.atlas = 'hcp_distal'; 
    end
end

switch op.atlas 

    %% 
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

        region_areas = {   'SMC',  {'1','2','3a','3b','4','6v','6d','43','55b','PEF','FEF','OP4','i6-8'};... % sensorimotor cortex... included operculum bc ecog strips can't get into operc
                        'vPMC', {'6r','FOP1'};... % ventral premotor... there are some elcs put into 'SMC' areas that are actually in vPMC... included operculum bc ecog strips can't get into operc
                        'STG', {'A4','A5','STGa','STV','TPOJ1' }; ... % superior temporal gyrus
                        'MFG',  {'8Av','8C','p9-46v'};... middle frontal gyrus... maybe also inf front sulcus
                        'IFG/IFS',  {'44','45','IFSp'};... % inferior frontal gyrus
                        'SMG/PF', {'PF','PFop'};... % supramarginal gyrus, operculum, ventral postcentral sulcus
                        'AG', {'PSL'};... % angular gyrus
        
                        'MTG', {'TE1a','TE1m','TE1p'};... middle temporal gyrus / TE
                        'STN', {'STN_associative_L','STN_motor_L','STN_motor_R' };...
                        'Thal', {'087_Thalamus_ventro_oralis_anterior_Voa_L','088_Thalamus_ventro_oralis_posterior_Vop_L','088_Thalamus_ventro_oralis_posterior_Vop_R',...
                                   '090_Thalamus_zentrolateralis_oralis_Zo_L','091_Thalamus_ventro_intermedius_internus_Vimi_R',...
                                   '094_Thalamus_ventro_intermedius_externus_Vime_L','094_Thalamus_ventro_intermedius_externus_Vime_R' };...
                        'GP', {'GPe_L','GPe_R','GPi_postparietal_R','GPi_premotor_R','GPi_sensorimotor_L','GPi_sensorimotor_L'};... % 
                        };

        op_out.atlas_var_names = {'HCPMMP1_label_1';'DISTAL_label_1'}; 
end

op_out.regiondef = table(region_areas(:,1), region_areas(:,2), 'Rownames', region_areas(:,1), 'VariableNames',...
                            {'region',          'areas'} );

op_out.nregions = height(op_out.regiondef); 

if exist('resp','var')
    resp_out = resp; 
    resp_out.region = cell(height(resp_out),1);
    resp_out = movevars(resp_out,'region','Before',op_out.atlas_var_names{1}); 
    for iregion = 1:op_out.nregions
        thisregion = op_out.regiondef.region{iregion};
            for iatlas = 1:length(op_out.atlas_var_names) 
                 atlas_var = op_out.atlas_var_names{iatlas}; 
                 elcs_in_this_region = ismember(resp_out{:,atlas_var},op_out.regiondef.areas{iregion}); 
                 resp_out.region(elcs_in_this_region) = {thisregion};
            end
    end
end


