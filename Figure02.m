%% Figure02 - Beta and FC fingerprints: shared backbone, distinct stable gradients
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions'));
fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('Figure02');
C=fig_colors();
fcdir=fullfile(cfg.result_dirs{12},'s011b_fc_vs_beta');
stabdir=fullfile(cfg.result_dirs{16},'s015a_beta_fc_stability');
gdir=fullfile(cfg.result_dirs{16},'s015b_gradient_stability');

G=readtable(fullfile(fcdir,'FC_vs_beta_group_summary.csv'),'Delimiter',',','TextType','string');
label_map = containers.Map( ...
    {'signed profile correlation','absolute profile correlation','17-network composition correlation','top ROI overlap'}, ...
    {'Signed profile','Absolute profile','17-network composition','Top-ROI overlap'});
labs = strings(height(G),1);
for i=1:height(G)
    key=char(G.group_metric(i));
    if isKey(label_map,key), labs(i)=label_map(key); else, labs(i)=strrep(string(key),'_',' '); end
end
fig=figure('Position',[100 100 390 270]); 
barh(1:height(G),G.group_profile_value,'FaceColor',[.45 .45 .45],'EdgeColor','none');
xlim([0 1]); yticks(1:height(G)); yticklabels(labs); xlabel('Group value');
title('Beta-FC similarity'); set(gca,'TickDir','out','YDir','reverse');
fig_export(fig,fullfile(paneldir,'panel_a_beta_fc_group_similarity.png'));

S=readtable(fullfile(fcdir,'FC_vs_beta_subject_summary.csv'),'Delimiter',',','TextType','string');
metrics={'mean_signed_profile_r','mean_absolute_profile_r','mean_network17_r','mean_topROI_overlap'};
fig=figure('Position',[100 100 390 280]); hold on;
X=[]; Y=[]; for k=1:numel(metrics), X=[X;k*ones(height(S),1)]; Y=[Y;S.(metrics{k})]; end
boxchart(X,Y,'BoxFaceColor',C.lightgray,'BoxEdgeColor',C.gray,'MarkerStyle','.'); ylim([0 1]); xticks(1:4); xticklabels({'Signed','Absolute','17-net','Top ROI'}); ylabel('Subject value'); title('Subject-level correspondence'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_b_subject_similarity_distributions.png'));

R=readtable(fullfile(stabdir,'beta_vs_FC_stability_subject_summary.csv'),'Delimiter',',','TextType','string');
R=R(R.split_name=="session",:);
fig_paired_cloud([R.mean_beta_stability R.mean_FC_stability],{'Beta','FC'}, ...
    'Raw fingerprint stability r','High-dimensional fingerprints', ...
    fullfile(paneldir,'panel_c_raw_fingerprint_stability.png'));

GS=readtable(fullfile(gdir,'beta_vs_FC_gradient_stability.csv'),'Delimiter',',','TextType','string');
G1=GS(GS.split_name=="session",:);
fig=figure('Position',[100 100 320 260]); cats=categorical(G1.modality); bar(cats,G1.gradient1_r,'FaceColor',[.45 .45 .45],'EdgeColor','none'); ylim([.85 1]); ylabel('Split/session G1 r'); title('Low-dimensional gradient stability'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_d_gradient1_stability.png'));

fig_brain_panel(fullfile(cfg.output_dir,'result_s023_beta_fc_gradient_templates','beta_G1_fixed_template.nii.gz'),fullfile(paneldir,'panel_e_beta_G1_template.png'),'Beta G1 fixed template','coolwarm',symmetric=true);
fig_brain_panel(fullfile(cfg.output_dir,'result_s023_beta_fc_gradient_templates','FC_G1_fixed_template.nii.gz'),fullfile(paneldir,'panel_f_FC_G1_template.png'),'FC G1 sign-aligned template','coolwarm',symmetric=true);
fig_brain_panel(fullfile(cfg.output_dir,'result_s023_beta_fc_gradient_templates','beta_minus_FC_G1_template.nii.gz'),fullfile(paneldir,'panel_g_beta_minus_FC_G1.png'),'Beta minus FC G1','coolwarm',symmetric=true);

fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure 2 layout guide\n\nRecommended layout: 2 rows. Row 1: a,b,c,d statistical panels. Row 2: e beta G1, f FC G1, g beta-FC G1 map. Keep brain maps same width and matched diverging colormap.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure 2 legend\n\nBeta and FC fingerprints share a strong backbone but diverge in low-dimensional organization. (a) Group beta-FC similarity metrics. (b) Subject-level distributions. (c) Raw 400-dimensional fingerprint stability shows FC exceeding beta. (d) In contrast, the first beta gradient is more reproducible than the first FC gradient. (e-g) Spatial maps of beta G1, FC G1, and their difference.\n");
fprintf('Figure02 complete: %s\n',figdir);
