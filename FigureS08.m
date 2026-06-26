%% FigureS08 - Spatial-control residual maps
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS08');
fig_brain_panel(fullfile(cfg.output_dir,'result_s032_spatial_confound_control','group_R2_residual_after_spatial_confounds.nii.gz'),fullfile(paneldir,'panel_a_controlled_R2_residual.png'),'Controlled R² residual','coolwarm',symmetric=true);
fig_brain_panel(fullfile(cfg.output_dir,'result_s032_spatial_confound_control','group_beta_minus_FC_G1_after_spatial_confounds.nii.gz'),fullfile(paneldir,'panel_b_controlled_beta_minus_FC_G1.png'),'Controlled beta-FC G1','coolwarm',symmetric=true);
fig_brain_panel(fullfile(cfg.output_dir,'result_s032_spatial_confound_control','group_zR2_after_zFC_and_spatial_confounds.nii.gz'),fullfile(paneldir,'panel_c_zR2_after_FC_and_confounds.png'),'zR2 after zFC+confounds','coolwarm',symmetric=true);
S=readtable(fullfile(cfg.output_dir,'result_s036_local_story_summary','local_story_extension_summary.csv'),'Delimiter',',','TextType','string');
labels={'R² residual','G1 diff','Joint'};
cats=categorical(labels);
cats=reordercats(cats,labels);
vals=[S.value(S.metric=="control_R2_residual_controlled_FWER_voxels") S.value(S.metric=="control_G1diff_controlled_FWER_voxels") S.value(S.metric=="control_joint_controlled_FWER_voxels")];
fig=figure('Position',[100 100 390 260]); bar(cats,vals,'FaceColor',[.35 .35 .35],'EdgeColor','none'); ylabel('FWER voxels'); title('Controlled significant voxels'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_d_controlled_FWER_extent.png'));
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S8 layout guide\n\nArrange three controlled residual maps plus FWER extent summary.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S8 legend\n\nSpatially controlled residual maps and significant voxel extents after nuisance adjustment.\n");
fprintf('FigureS08 complete\n');
