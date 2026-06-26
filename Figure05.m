%% Figure05 - Continuous G1 axis: cortical sources and tract anatomy
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions'));
fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('Figure05');
C=fig_colors();

fig_brain_panel(fullfile(cfg.output_dir,'result_s024_beta_fc_gradient_subject_inference','subject_mean_beta_minus_FC_G1.nii.gz'),fullfile(paneldir,'panel_a_beta_minus_FC_G1_axis.png'),'Beta-FC G1 axis','coolwarm',symmetric=true);
B=readtable(fullfile(cfg.output_dir,'result_s034_g1_binned_cortical_contribution','G1_bin_summary.csv'),'Delimiter',',','TextType','string');

fig=figure('Position',[100 100 390 280]); hold on;
plot(B.bin,B.mean_zR2,'-o','Color',C.blue,'MarkerFaceColor','w','MarkerSize',3.5);
plot(B.bin,B.mean_zFC_RMS,'-o','Color',C.gray,'MarkerFaceColor','w','MarkerSize',3.5);
plot(B.bin,B.mean_R2_residual,'-o','Color',C.red,'MarkerFaceColor','w','MarkerSize',3.5);
xlabel('Beta-FC G1 bin'); ylabel('Mean z-score / residual'); legend({'zR²','zFC RMS','R² residual'},'Location','southoutside','Orientation','horizontal'); title('Prediction structure along G1'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_b_prediction_along_G1.png'));

fig=figure('Position',[100 100 360 270]); hold on;
plot(B.bin,B.mean_R2_residual,'-o','Color',C.red,'MarkerFaceColor','w','MarkerSize',3.5);
plot(B.bin,B.mean_R2_residual_controlled,'-o','Color',C.teal,'MarkerFaceColor','w','MarkerSize',3.5);
plot(B.bin,B.mean_zR2_after_zFC_spatial_confounds,'-o','Color',C.blue,'MarkerFaceColor','w','MarkerSize',3.5);
xlabel('Beta-FC G1 bin'); ylabel('Mean residual'); legend({'FC residual','Controlled','zR² | zFC + conf'},'Location','southoutside','Orientation','horizontal'); title('Residual axis survives control'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_c_controlled_residual_along_G1.png'));

N=readtable(fullfile(cfg.output_dir,'result_s034_g1_binned_cortical_contribution','G1_bin_network_contribution.csv'),'Delimiter',',','TextType','string');
nets=unique(N.network_name,'stable'); M=zeros(numel(nets),max(N.bin));
for i=1:height(N), M(nets==N.network_name(i),N.bin(i))=N.fraction_abs_DeltaR2(i); end
fig=figure('Position',[100 100 470 360]); imagesc(M); colormap(fig_sequential_colormap(256,'bluegreen')); cb=colorbar; cb.Label.String='Fraction |ΔR²|'; yticks(1:numel(nets)); yticklabels(nets); xticks(1:max(N.bin)); xlabel('Beta-FC G1 bin'); title('Cortical network contribution');
set(gca,'TickDir','out','FontSize',6.5); fig_export(fig,fullfile(paneldir,'panel_d_network_contribution_heatmap.png'));

J=readtable(fullfile(cfg.output_dir,'result_s034_g1_binned_cortical_contribution','G1_bin_JHU_composition.csv'),'Delimiter',',','TextType','string');
tot=groupsummary(J,'label_name','sum','direct_voxels'); [~,ord]=sort(tot.sum_direct_voxels,'descend'); keep=tot.label_name(ord(1:min(10,height(tot))));
JM=zeros(numel(keep),max(J.bin));
for i=1:height(J)
    k=find(keep==J.label_name(i),1); if ~isempty(k), JM(k,J.bin(i))=J.fraction_within_bin(i); end
end
keep_labels=regexprep(string(keep),'\s*\([^)]*\)','');
keep_labels=strrep(keep_labels,'Posterior thalamic radiation R','Posterior thalamic rad. R');
fig=figure('Position',[100 100 520 330]); imagesc(JM); colormap(flipud(gray(256))); cb=colorbar; cb.Label.String='Fraction within bin'; yticks(1:numel(keep)); yticklabels(keep_labels); xticks(1:max(J.bin)); xlabel('Beta-FC G1 bin'); title('JHU composition across G1');
set(gca,'TickDir','out','FontSize',6.5); fig_export(fig,fullfile(paneldir,'panel_e_JHU_composition_heatmap.png'));

H=readtable(fullfile(cfg.output_dir,'result_s035_hemispheric_symmetry','hemispheric_symmetry_map_correlations.csv'),'Delimiter',',','TextType','string');
map_labels=strrep(string(H.map_name),'R2','R²');
cats=categorical(map_labels);
cats=reordercats(cats,cellstr(map_labels));
fig=figure('Position',[100 100 360 270]); bar(cats,H.LR_correlation,'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylim([0 1]); xtickangle(35); ylabel('Left-right mirror r'); title('Bilateral organization'); set(gca,'TickDir','out');
fig_export(fig,fullfile(paneldir,'panel_f_hemispheric_symmetry.png'));

fig=figure('Position',[100 100 620 280]); ax=axes(fig,'Position',[.035 .08 .93 .84]); axis(ax,'off'); hold(ax,'on');
rectangle(ax,'Position',[.03 .28 .34 .50],'Curvature',.05,'FaceColor',[.96 .97 .98],'EdgeColor',C.gray,'LineWidth',.7);
rectangle(ax,'Position',[.63 .28 .34 .50],'Curvature',.05,'FaceColor',[.98 .96 .94],'EdgeColor',C.gray,'LineWidth',.7);
text(ax,.20,.70,'Marginal FC backbone','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
text(ax,.20,.52,sprintf('Shared beta-FC organization\nhigh signed profile correlation'), ...
    'FontSize',7,'HorizontalAlignment','center','Color',[.20 .20 .20]);
text(ax,.80,.70,'Conditional encoding axis','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
text(ax,.80,.52,sprintf('Low marginal FC\nhigh FC-controlled R^2 residual'), ...
    'FontSize',7,'HorizontalAlignment','center','Color',[.20 .20 .20]);
annotation(fig,'arrow',[.42 .58],[.56 .56],'LineWidth',.9,'HeadLength',7,'HeadWidth',7,'Color',C.gray);
text(ax,.20,.22,'Callosal / capsule end','FontSize',7,'HorizontalAlignment','center','Color',C.gray);
text(ax,.80,.22,'Posterior thalamic / optic radiation end','FontSize',7,'HorizontalAlignment','center','Color',C.gray);
text(ax,.50,.33,'G1','FontSize',7,'FontWeight','bold','HorizontalAlignment','center','Color',C.gray);
fig_export(fig,fullfile(paneldir,'panel_g_conceptual_model.png'));

fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure 5 layout guide\n\nRecommended layout: 2 rows. Row 1: a G1 map, b prediction curves, c controlled residual curves. Row 2: d network heatmap, e JHU composition heatmap, f hemispheric symmetry, g conceptual model. Heatmaps may need extra width because tract/network names are long.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure 5 legend\n\nA continuous conditional encoding axis organizes white-matter BOLD. (a) Beta minus FC G1 axis. (b) Along the axis, high-G1 bins show low marginal FC but high FC-controlled prediction residuals. (c) This axis remains after nuisance controls. (d) Cortical network contribution shifts along the axis, with visual/default contributions increasing toward the high-G1 end. (e) JHU tract composition shifts from callosal/internal-capsule structures toward posterior thalamic/optic radiation and corona-radiata structures. (f) Key maps show bilateral mirror symmetry. (g) Conceptual distinction between marginal FC backbone and conditional cortical encoding.\n");
fprintf('Figure05 complete: %s\n',figdir);
