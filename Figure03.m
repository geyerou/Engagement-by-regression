%% Figure03 - FC-controlled local prediction converges with beta-FC G1
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions'));
fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('Figure03');
COL=fig_colors();

fig_brain_panel(fullfile(cfg.output_dir,'result_s021_local_r2_fc','group_mean_local_zR2.nii.gz'),fullfile(paneldir,'panel_a_local_zR2.png'),'Local zR²','viridis');
fig_brain_panel(fullfile(cfg.output_dir,'result_s021_local_r2_fc','group_mean_local_zFC_RMS.nii.gz'),fullfile(paneldir,'panel_b_local_zFC_RMS.png'),'Local zFC RMS','viridis');
fig_brain_panel(fullfile(cfg.output_dir,'result_s021_local_r2_fc','group_mean_zR2_residualized_local_zFC.nii.gz'),fullfile(paneldir,'panel_c_R2_residual_after_FC.png'),'R² residual after FC','coolwarm',symmetric=true);

C=readtable(fullfile(cfg.output_dir,'result_s022_local_r2_fc_inference','local_R2_residual_clusters.csv'),'Delimiter',',','TextType','string');
sig=C(C.permutation_p<.05,:);
fig=figure('Position',[100 100 340 270]); 
bar(categorical(["Positive","Negative"]),[sum(sig.direction=="positive") sum(sig.direction=="negative")],'FaceColor',COL.gray,'EdgeColor','none');
ylabel('FWER-significant clusters'); title('Signed local R² residual clusters'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_d_significant_cluster_counts.png'));

fig_brain_panel(fullfile(cfg.output_dir,'result_s024_beta_fc_gradient_subject_inference','subject_mean_beta_minus_FC_G1.nii.gz'),fullfile(paneldir,'panel_e_subject_mean_beta_minus_FC_G1.png'),'Subject mean beta-FC G1','coolwarm',symmetric=true);

L=load(fullfile(cfg.output_dir,'result_s021_local_r2_fc','local_r2_fc_group_data.mat'),'group_mean_residual');
G=load(fullfile(cfg.output_dir,'result_s024_beta_fc_gradient_subject_inference','gradient_subject_inference.mat'),'mean_difference');
r=corr(double(L.group_mean_residual),double(G.mean_difference));
fig_density_scatter(L.group_mean_residual,G.mean_difference, ...
    'Local R² residual','Beta minus FC G1',sprintf('r = %.3f',r), ...
    fullfile(paneldir,'panel_f_raw_residual_vs_G1_scatter.png'));

K=load(fullfile(cfg.output_dir,'result_s032_spatial_confound_control','spatial_confound_control.mat'),'group_R2_residual_controlled','group_G1diff_controlled','r2_pfwer','g1_pfwer','C');
rc=corr(double(K.group_R2_residual_controlled),double(K.group_G1diff_controlled));
fig_density_scatter(K.group_R2_residual_controlled,K.group_G1diff_controlled, ...
    'Controlled R² residual','Controlled beta-FC G1',sprintf('r = %.3f',rc), ...
    fullfile(paneldir,'panel_g_controlled_scatter.png'));

joint=(K.r2_pfwer(:)<.05) & (K.g1_pfwer(:)<.05) & sign(K.group_R2_residual_controlled)==sign(K.group_G1diff_controlled);
write_map_like(single(joint),K.C.wm_idx,cfg.wm_mask_final,fullfile(paneldir,'panel_h_joint_same_sign_controlled_map.nii.gz'));
fig_brain_panel(fullfile(paneldir,'panel_h_joint_same_sign_controlled_map.nii.gz'),fullfile(paneldir,'panel_h_joint_same_sign_controlled_map.png'),'Joint controlled voxels','Reds',kind='roi');

fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure 3 layout guide\n\nRecommended layout: 2 rows x 4 columns. Row 1: a local zR2, b local zFC, c FC-controlled R2 residual, d signed cluster counts. Row 2: e beta-FC G1, f raw residual-G1 scatter, g controlled scatter, h joint same-sign controlled map.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure 3 legend\n\nFC-controlled local prediction residuals converge with the beta-FC gradient difference. (a-c) Group maps of local prediction, local marginal FC strength, and R2 residual after controlling local FC. (d) Number of signed cluster-mass FWER significant residual clusters. (e) Mean beta minus FC G1. (f) Spatial correspondence between local R2 residual and beta-FC G1. (g) The correspondence remains after anatomical/signal/spatial nuisance control. (h) Same-sign voxels significant for both controlled effects.\n");
fprintf('Figure03 complete: %s\n',figdir);
