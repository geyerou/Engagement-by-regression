%% Figure04 - Spatial confounds and spatial-autocorrelation null
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions'));
fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('Figure04');

fig_brain_panel(fullfile(cfg.output_dir,'result_s031_spatial_confound_maps','distance_to_group_GM_mm.nii.gz'),fullfile(paneldir,'panel_a_distance_to_GM.png'),'Distance to GM','magma');
fig_brain_panel(fullfile(cfg.output_dir,'result_s031_spatial_confound_maps','distance_to_WM_boundary_mm.nii.gz'),fullfile(paneldir,'panel_b_distance_to_WM_boundary.png'),'Distance to WM boundary','magma');
fig_brain_panel(fullfile(cfg.output_dir,'result_s031_spatial_confound_maps','subject_WM_mask_prevalence.nii.gz'),fullfile(paneldir,'panel_c_WM_prevalence.png'),'WM prevalence','viridis',vmax=1);
fig_brain_panel(fullfile(cfg.output_dir,'result_s031_spatial_confound_maps','mean_extracted_WM_temporal_SD.nii.gz'),fullfile(paneldir,'panel_d_temporal_SD.png'),'Temporal SD','viridis');

S=readtable(fullfile(cfg.output_dir,'result_s032_spatial_confound_control','spatial_confound_control_summary.csv'),'Delimiter',',','TextType','string');
metrics=["map_r_R2resid_G1diff_before_control","map_r_R2resid_G1diff_after_control"];
vals=arrayfun(@(m)S.value(S.metric==m),metrics);
fig=figure('Position',[100 100 300 260]); bar(categorical(["Before","After"]),vals,'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylim([0 .6]); ylabel('Spatial r'); title('Effect survives nuisance control'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_e_control_before_after_r.png'));

N=load(fullfile(cfg.output_dir,'result_s033_spatial_autocorr_null','spatial_block_null.mat'));
fig=figure('Position',[100 100 330 260]); histogram(N.null_ctl_r,40,'FaceColor',[.70 .70 .70],'EdgeColor','none'); hold on; xline(N.obs_ctl_r,'Color',[.85 .25 .20],'LineWidth',1.2);
xlabel('Block-null spatial r'); ylabel('Permutations'); title(sprintf('Controlled r block p=%.4f',N.p_ctl_r)); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_f_spatial_block_null_r.png'));

fig=figure('Position',[100 100 330 260]); histogram(N.null_same_sign,40,'FaceColor',[.70 .70 .70],'EdgeColor','none'); hold on; xline(N.obs_same_sign,'Color',[.85 .25 .20],'LineWidth',1.2);
xlabel('Same-sign joint voxels'); ylabel('Permutations'); title(sprintf('Joint block p=%.4f',N.p_same_sign)); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_g_spatial_block_null_joint.png'));

B=readtable(fullfile(cfg.output_dir,'result_s034_g1_binned_cortical_contribution','G1_bin_summary.csv'),'Delimiter',',','TextType','string');
fig=figure('Position',[100 100 380 270]);
yyaxis left;
plot(B.bin,B.mean_dist_GM_mm,'-o','Color',[.12 .31 .54],'MarkerFaceColor','w','MarkerSize',3.5);
ylabel('Distance to GM (mm)');
ax=gca; ax.YColor=[.12 .31 .54];
yyaxis right;
plot(B.bin,B.mean_WM_prevalence,'-o','Color',[0 .48 .42],'MarkerFaceColor','w','MarkerSize',3.5);
ylabel('WM prevalence');
ax=gca; ax.YAxis(2).Color=[0 .48 .42];
xlabel('Beta-FC G1 bin'); title('High G1 is not GM-proximal'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_h_G1_bin_contamination_diagnostics.png'));

fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure 4 layout guide\n\nRecommended layout: 2 rows x 4 columns. Row 1: a-d nuisance/proximity maps. Row 2: e control before-after, f block-null map correlation, g block-null joint overlap, h G1-bin contamination diagnostics.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure 4 legend\n\nSpatial control analyses argue against gray-matter contamination or spatial autocorrelation as sufficient explanations. (a-d) Nuisance maps for GM distance, WM boundary distance, subject WM prevalence, and temporal signal SD. (e) R2 residual to beta-FC G1 spatial correspondence before and after nuisance control. (f-g) Volumetric WM spatial block-permutation null distributions. (h) G1-bin diagnostics show that the high-G1/high-residual end is not simply closer to GM and has high WM prevalence.\n");
fprintf('Figure04 complete: %s\n',figdir);
