%% Figure01 - Strict WM ridge prediction and circular-shift null
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions'));
fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('Figure01');
C=fig_colors();

% a. Encoding model schematic
fig_ridge_encoding_schematic(cfg.schaefer400_file,cfg.wm_mask_final, ...
    fullfile(paneldir,'panel_a_pipeline.png'));

% b-c. Brain maps
fig_brain_panel(cfg.wm_mask_final,fullfile(paneldir,'panel_b_strict_wm95_mask.png'),'Strict WM95 mask','Reds',kind='roi');
fig_brain_panel(fullfile(cfg.result_dirs{17},'group_mean_R2_map.nii.gz'),fullfile(paneldir,'panel_c_group_mean_R2.png'),'Group mean R^2','viridis',vmax=.08);

% d. Subject mean R2 distribution
M=readtable(fullfile(cfg.result_dirs{19},'model_comparison_table.csv'),'Delimiter',',','TextType','string');
Z=M(M.model=="ridge_400_zero_lag",:);
fig=figure('Position',[100 100 320 260]); histogram(Z.mean_R2,18,'FaceColor',[.55 .55 .55],'EdgeColor','none'); xline(mean(Z.mean_R2),'Color',C.red,'LineWidth',1.2);
xlabel('Subject mean R²'); ylabel('Subjects'); title('Zero-lag GM-to-WM prediction'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_d_subject_mean_R2_distribution.png'));

% e. 400 ROI versus 17-network model
models=["ridge_400_zero_lag","ridge_network17"]; vals=zeros(height(unique(Z.subject_id)),2);
S=unique(M.subject_id,'stable');
for i=1:numel(S)
    for j=1:2, vals(i,j)=M.mean_R2(M.subject_id==S(i) & M.model==models(j)); end
end
fig_paired_cloud(vals,{'400 ROI','17 networks'},'Subject mean R²', ...
    'Fine-grained GM predictors matter',fullfile(paneldir,'panel_e_400roi_vs_17network.png'));

% f. Full s014 circular-shift null extent
qfiles=dir(fullfile(cfg.result_dirs{15},'*_null_qmap_R2.nii.gz'));
wm=niftiread(cfg.wm_mask_final)>0; qfrac=zeros(numel(qfiles),1); pfrac=zeros(numel(qfiles),1);
for i=1:numel(qfiles)
    q=niftiread(fullfile(qfiles(i).folder,qfiles(i).name));
    sid=extractBetween(qfiles(i).name,'sub-','_null_qmap_R2');
    p=niftiread(fullfile(qfiles(i).folder,sprintf('sub-%s_null_pmap_R2.nii.gz',sid{1})));
    qfrac(i)=mean(q(wm)<.05); pfrac(i)=mean(p(wm)<.05);
end
fig=figure('Position',[100 100 330 270]); hold on;
histogram(qfrac,18,'FaceColor',[.55 .55 .55],'EdgeColor','none'); xline(mean(qfrac),'Color',C.red,'LineWidth',1.2);
xlabel('Fraction WM voxels q_null < 0.05'); ylabel('Subjects'); title('Circular-shift null, n=81'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_f_s014_null_extent.png'));

layout=sprintf(['# Figure 1 layout guide\n\n' ...
'Recommended layout: 2 rows. Row 1: a wide pipeline panel spanning two columns, then b strict WM mask and c group mean R2. Row 2: d subject distribution, e model comparison, f circular-shift null extent.\n\n' ...
'Panel files are in `panels/`. Use 9 pt Arial labels; panel letters 10-11 pt bold.\n']);
legend=sprintf(['# Figure 1 legend\n\n' ...
'Strict white-matter cortical prediction model. (a) Analysis workflow from HCP resting-state fMRI through Schaefer-400 gray-matter predictors, voxelwise ridge prediction in strict WM95, and fingerprint/gradient analyses. (b) Strict group WM95 mask used for all primary analyses. (c) Group mean cross-validated R2 map for the zero-lag 400-ROI ridge model. (d) Subject distribution of mean WM prediction accuracy. (e) Paired comparison of fine-grained Schaefer-400 predictors versus 17-network predictors. (f) Full circular-shift null validation from s014, showing the fraction of WM voxels surviving voxelwise FDR q<0.05 per subject.\n']);
fig_write_text(fullfile(figdir,'layout_guide.md'),layout);
fig_write_text(fullfile(figdir,'legend.md'),legend);
fprintf('Figure01 complete: %s\n',figdir);
