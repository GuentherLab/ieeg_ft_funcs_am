%%%% plot electrodes meeting certain criteria on brain 


function [op] = brainplot_top_electrodes(resp,op)

field_default('op','tuning_param',{})              % if tuning_param is empty, plot all electrodes without sorting
field_default('op','inclusion_mode','threshold');   % options: 'threshold','proportion'
    field_default('op','include_proportion')        % only has effect if inclusion_mode == proportion
    field_default('op','include_thresh',0.05)     % only has effect if inclusion_mode == threshold
field_default('op','struct_to_plot','ctx') % options: ctx, stn, thal
field_default('op','snap_to_surf',true) % cortex only - if true, project eletrodes to nearest point on ctx surface
    field_default('op','snap_offset_x',-1) % .... if snapping, offset of -1 should be enough to have points entirely above ctx surface (in L hem)
field_default('op','hemi','L') % only applies to subcortical..... options: L, R
field_default('op','view_angle',  [-90, 0] ) % use [-90, 0] for straight-on lateral left hemisphere
field_default('op','elc_types_to_plot',{'ECOG'});

field_default('op','marker_size',40)
field_default('op','plotcolor','r')
field_default('op','also_plot_nonsgnf_elcs',true) % if true, plot elcs not meeting criteria in different color and size
    field_default('op','marker_size_nonsgn',5) 
    field_default('op','plotcolor_nonsgn',[0.3 0.3 0.3])



% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')





%% Configuration Variables and Paths
% PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures';
% % % % % % % % % PATH_DATA='/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures/data';
field_default('op','path_average_mni','Z:\DBS\DBS_subject_lists/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexLowRes_15000V.mat')
field_default('op','path_subcort_atlas','/Volumes/Nexus/Resources/STN-Atlas/atlas_index.mat')
field_default('op','path_subcort_atlas_vim','/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/atlas_index.mat')


% cd(PATH_ANALYSIS)
% electrode = readtable('data/A01_DBS_aper_coord_dx.tsv','Delimiter', '\t', 'TreatAsEmpty', 'NA','FileType','text');

%loading cortical reconstructions
average_mni = load(op.path_average_mni);

% subcort = load(op.path_subcort_atlas);
% subcort_vim = load(op.path_subcort_atlas_vim);
% nii_vimi = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/91.nii.gz');
% nii_vime = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/94.nii.gz');
% nii_vimip = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/104.nii.gz');
% nii_vimep = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/122.nii.gz');
% 
% subcort_vimi_lh_fv = ea_nii2fv(nii_vimi);
% subcort_vime_lh_fv = ea_nii2fv(nii_vime);
% subcort_vimip_lh_fv = ea_nii2fv(nii_vimip);
% subcort_vimep_lh_fv = ea_nii2fv(nii_vimep);

% % % %loading VL posterior ventral from Morel atlas
% nii_vlpv = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/MorelAtlasICBM2009b (Jakab 2008)/lh/VLpv.nii.gz');
% subcort_vlpv_lh_fv = ea_nii2fv(nii_vlpv);

color_et_ecog = '#C4604F';% #ET ECoG
color_pd_ecog = '#6F67A6';% #PD ECoG
color_ep_seeg = '#8A4F80';% #EP sEEG
color_pd_stn = '#F7924A';% #PD STN
color_pd_gpi = '#F9BD00';% #PD GPi
color_et_vim = '#36A5D1';% #ET VIM
color_ep_cm = '#9EB859';% #EP CM


n_elc = height(resp);
if isempty(op.tuning_param)
    sgn_rows = true(height(resp),1); % include all electrodes if not inclusion variable was specified
else
    switch op.inclusion_mode
        case 'thresh'
            sgn_rows = resp{:,op.tuning_param} < op.include_thresh;
        case 'proportion'
            [~, rows_ranked] = sort(resp{:,op.tuning_param});
            sgn_rows = rows_ranked( 1:round(op.include_proportion * n_elc) ); 
        otherwise
            error('unknown inclusion mode')
    end
end

%% load brain surfaces
switch op.struct_to_plot
    case 'ctx'
        %loading cortical reconstructions
        average_mni = load(op.path_average_mni);
    case 'stn'
        subcort_stn = load(PATH_STN_ATLAS);
    case 'thal'
        subcort_vim = load(op.path_subcort_atlas_vim);
        % nii_vimi = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/91.nii.gz');
        % nii_vime = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/94.nii.gz');
        % nii_vimip = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/104.nii.gz');
        % nii_vimep = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/122.nii.gz');
        
        % 
        % subcort_vimi_lh_fv = ea_nii2fv(nii_vimi);
        % subcort_vime_lh_fv = ea_nii2fv(nii_vime);
        % subcort_vimip_lh_fv = ea_nii2fv(nii_vimip);
        % subcort_vimep_lh_fv = ea_nii2fv(nii_vimep);
        
        % % % %loading VL posterior ventral from Morel atlas
        % nii_vlpv = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/MorelAtlasICBM2009b (Jakab 2008)/lh/VLpv.nii.gz');
        % subcort_vlpv_lh_fv = ea_nii2fv(nii_vlpv);
end


%% make brainplot

switch op.hemi
    case 'L'
        op.hemi_number = 2;
    case 'R'
        op.hemi_number = 1; 
end


% close all

% % % % % % % % % % % % % % % hfig = figure('WindowState','maximized');
% % % % % % % % % % % % % % % patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
% % % % % % % % % % % % % % % 'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
% % % % % % % % % % % % % % % 'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
% % % % % % % % % % % % % % % hold on




% rows_to_plot & any(strcmp(resp.type,op.elc_types_to_plot), 2);



% % % % if op.also_plot_nonsgnf_elcs
% % % %     elc_to_plot = resp(~rows_to_plot,{'mni_x','mni_y','mni_z'}); 
% % % % 
% % % %     xyz_to_plot_nonsnapped = [elc_to_plot.mni_x, elc_to_plot.mni_y, elc_to_plot.mni_z];
% % % %     if op.snap_to_surf
% % % %         [~, surfpoint_idx] = min(pdist2(xyz_to_plot_nonsnapped,average_mni.Vertices), [], 2); % find nearest surf points
% % % %         xyz_to_plot = average_mni.Vertices(surfpoint_idx,:); 
% % % %     elseif ~op.snap_to_surf
% % % %         xyz_to_plot = xyz_to_plot_nonsnapped;
% % % %     end
% % % % 
% % % %     hscat = scatter3(xyz_to_plot(:,1) + op.snap_offset_x, xyz_to_plot(:,2), xyz_to_plot(:,3), 'filled',...
% % % %   'MarkerFaceAlpha',1,'MarkerFaceColor',op.plotcolor_nonsgnf,'MarkerEdgeColor','k','LineWidth',0.01);
% % % %     hscat.SizeData = 20;
% % % %     % % % % % % % % % % % scalebar(0,70,-50, 10, 'mm')
% % % % end
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % elc_to_plot = resp(rows_to_plot,{'mni_x','mni_y','mni_z'}); 

%%

% xyz_to_plot_nonsnapped = [elc_to_plot.mni_x, elc_to_plot.mni_y, elc_to_plot.mni_z];
% if op.snap_to_surf
%     [~, surfpoint_idx] = min(pdist2(xyz_to_plot_nonsnapped,average_mni.Vertices), [], 2); % find nearest surf points
%     xyz_to_plot = average_mni.Vertices(surfpoint_idx,:); 
% elseif ~op.snap_to_surf
%     xyz_to_plot = xyz_to_plot_nonsnapped;
% end
% 
% hscat = scatter3(xyz_to_plot(:,1) + op.snap_offset_x, xyz_to_plot(:,2), xyz_to_plot(:,3), 'filled',...
%   'MarkerFaceAlpha',1,'MarkerFaceColor',op.plotcolor,'MarkerEdgeColor','k','LineWidth',0.01);
% hscat.SizeData = 60;
% set(gcf, 'Color', [1 1 1]); % white backgroud
% view(-90,0)
% axis off; axis equal
% camlight('headlight','infinite');
% % % % % % % % % % % % scalebar(0,70,-50, 10, 'mm')




switch op.struct_to_plot
    case 'ctx'
        rows_to_plot_sgn = sgn_rows & string(resp.type)=="ECOG"; 
        elc_to_plot = resp(rows_to_plot_sgn,{'mni_x','mni_y','mni_z'}); 
        xyz_to_plot_nonsnapped = [elc_to_plot.mni_x, elc_to_plot.mni_y, elc_to_plot.mni_z];

        rows_to_plot_nonsgn = ~sgn_rows & strcmp((resp.type), op.elc_types_to_plot); 
        elc_to_plot_nonsgn = resp(rows_to_plot_nonsgn,{'mni_x','mni_y','mni_z'}); 
        xyz_to_plot_nonsnapped_nonsgn = [elc_to_plot_nonsgn.mni_x, elc_to_plot_nonsgn.mni_y, elc_to_plot_nonsgn.mni_z];

        % shift electrodes so that they aren't covered by the brain surface
        %%% gets applied after snapping to surface
        if op.snap_to_surf
            [~, surfpoint_idx] = min(pdist2(xyz_to_plot_nonsnapped,average_mni.Vertices), [], 2); % find nearest surf points
            xyz_to_plot = average_mni.Vertices(surfpoint_idx,:); 
            [~, surfpoint_idx] = min(pdist2(xyz_to_plot_nonsnapped_nonsgn,average_mni.Vertices), [], 2); % find nearest surf points
            xyz_to_plot_nonsgn = average_mni.Vertices(surfpoint_idx,:); 
        elseif ~op.snap_to_surf
            xyz_to_plot = xyz_to_plot_nonsnapped;
            xyz_to_plot_nonsgn = xyz_to_plot_nonsnapped_nonsgn;
        end

        hpatch = patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
            'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
            'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);


    case 'stn'
        rows_to_plot_sgn = sgn_rows & contains(resp.type,{'DBS';'MACRO'}) & contains(resp.DISTAL_label_1,{'STN_'}) & contains(resp.DISTAL_label_1,{['_',op.hemi]});
        elc_to_plot = resp(rows_to_plot_sgn,{'mni_x','mni_y','mni_z'}); 
        xyz_to_plot = [elc_to_plot.mni_x, elc_to_plot.mni_y, elc_to_plot.mni_z];

        rows_to_plot_nonsgn = ~sgn_rows & contains(resp.type,{'DBS';'MACRO'}) & contains(resp.DISTAL_label_1,{'STN_'}) & contains(resp.DISTAL_label_1,{['_',op.hemi]});
        elc_to_plot_nonsgn = resp(rows_to_plot_nonsgn,{'mni_x','mni_y','mni_z'}); 
        xyz_to_plot_nonsgn = [elc_to_plot_nonsgn.mni_x, elc_to_plot_nonsgn.mni_y, elc_to_plot_nonsgn.mni_z];

        hpatch = patch('vertices', subcort_stn.atlases.fv{1,op.hemi_number}.vertices, 'faces', subcort_stn.atlases.fv{1,op.hemi_number}.faces,...
            'FaceColor', [.7 .6 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
            'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);

    case 'thal' %%%%%% need to find an appropriate VIM atlas before using this option
        rows_to_plot_sgn = sgn_rows & contains(resp.type,{'DBS';'MACRO'});
        rows_to_plot_sgn = rows_to_plot_sgn & contains(resp.MOREL_label_1,{'Thalamus'}) & contains(resp.DISTAL_label_1,{['_',op.hemi]});
        elc_to_plot = resp(rows_to_plot_sgn,{'mni_x','mni_y','mni_z'}); 
        xyz_to_plot = [elc_to_plot.mni_x, elc_to_plot.mni_y, elc_to_plot.mni_z];

        % hpatch = patch('vertices', subcort_vim.atlases.fv{???????,op.hemi_number}.vertices, 'faces', subcort_vim.atlases.fv{???????,op.hemi_number}.faces,...
        %     'FaceColor', [.7 .6 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
        %     'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);   


    case 'gp'  %%%%%% need to find an appropriate VIM atlas before using this option
        rows_to_plot_sgn = sgn_rows & contains(resp.type,{'DBS';'MACRO'}) & contains(resp.DISTAL_label_1,{'???'}) & contains(resp.DISTAL_label_1,{['_',op.hemi]});
        elc_to_plot = resp(rows_to_plot_sgn,{'mni_x','mni_y','mni_z'}); 
        xyz_to_plot = [elc_to_plot.mni_x, elc_to_plot.mni_y, elc_to_plot.mni_z];

        rows_to_plot_nonsgn = ~sgn_rows & contains(resp.type,{'DBS';'MACRO'}) & contains(resp.DISTAL_label_1,{'???'}) & contains(resp.DISTAL_label_1,{['_',op.hemi]});
        elc_to_plot_nonsgn = resp(rows_to_plot_nonsgn,{'mni_x','mni_y','mni_z'}); 
        xyz_to_plot_nonsgn = [elc_to_plot_nonsgn.mni_x, elc_to_plot_nonsgn.mni_y, elc_to_plot_nonsgn.mni_z];

        % hpatch = patch('vertices', subcort_vim.atlases.fv{???????,op.hemi_number}.vertices, 'faces', subcort_vim.atlases.fv{???????,op.hemi_number}.faces,...
        %     'FaceColor', [.7 .6 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
        %     'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);   

end

hold on

if  op.also_plot_nonsgnf_elcs
    hscat_non_sgnf = scatter3(xyz_to_plot_nonsgn(:,1) + op.snap_offset_x, xyz_to_plot_nonsgn(:,2), xyz_to_plot_nonsgn(:,3), 'filled',...
       'MarkerFaceAlpha',1,'MarkerFaceColor',op.plotcolor_nonsgn,'MarkerEdgeColor','k','LineWidth',0.01);
    hscat_non_sgnf.SizeData = op.marker_size_nonsgn;
end

hscat_sgnf = scatter3(xyz_to_plot(:,1) + op.snap_offset_x, xyz_to_plot(:,2), xyz_to_plot(:,3), 'filled',...
   'MarkerFaceAlpha',1,'MarkerFaceColor',op.plotcolor,'MarkerEdgeColor','k','LineWidth',0.01);
hscat_sgnf.SizeData = op.marker_size;

set(gcf, 'Color', [1 1 1]); % white backgroud
view(op.view_angle(1),op.view_angle(2))
axis off; axis equal
camlight('headlight','infinite');


titlestr = op.tuning_param; 
title(titlestr,'interpreter', 'none')

% print(gcf,[PATH_ANALYSIS 'qqq.png'],'-dpng','-r300')

