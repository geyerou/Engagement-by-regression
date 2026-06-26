%% FigureS06 - JHU localization
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS06');
J=readtable(fullfile(cfg.output_dir,'result_s025_jhu_localization','JHU_local_metric_summary.csv'),'Delimiter',',','TextType','string');
J=J(J.direct_WM_voxels>=20,:);
[~,op]=sort(J.direct_mean_R2_residual,'descend'); P=J(op(1:min(12,height(J))),:);
P=P(1:min(8,height(P)),:);
fig=figure('Position',[100 100 620 330]); barh(categorical(local_short_labels(P.label_name)),P.direct_mean_R2_residual,'FaceColor',[.85 .25 .20],'EdgeColor','none'); set(gca,'YDir','reverse','TickDir','out','FontSize',8); xlabel('Mean R² residual'); title('Top positive JHU tracts'); fig_export(fig,fullfile(paneldir,'panel_a_top_positive_R2_residual_tracts.png'));
[~,on]=sort(J.direct_mean_R2_residual,'ascend'); N=J(on(1:min(12,height(J))),:);
N=N(1:min(8,height(N)),:);
fig=figure('Position',[100 100 620 330]); barh(categorical(local_short_labels(N.label_name)),N.direct_mean_R2_residual,'FaceColor',[.33 .44 .68],'EdgeColor','none'); set(gca,'YDir','reverse','TickDir','out','FontSize',8); xlabel('Mean R² residual'); title('Top negative JHU tracts'); fig_export(fig,fullfile(paneldir,'panel_b_top_negative_R2_residual_tracts.png'));
fig=figure('Position',[100 100 330 280]); scatter(J.direct_mean_R2_residual,J.direct_mean_beta_minus_FC_G1,40,J.direct_WM_voxels,'filled'); colorbar; xlabel('Mean R² residual'); ylabel('Mean beta-FC G1'); title('Tract-level convergence'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_c_JHU_R2residual_vs_G1.png'));
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S6 layout guide\n\nArrange positive tracts, negative tracts, and tract-level R2 residual versus beta-FC G1 scatter.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S6 legend\n\nDirect-overlap JHU localization of local R2 residual and beta-FC G1 effects within strict WM95.\n");
fprintf('FigureS06 complete\n');

function labels=local_short_labels(labels)
labels=string(labels);
labels=strrep(labels,' (include optic radiation)','');
labels=regexprep(labels,'Sagittal stratum.* ([LR])$','Sagittal stratum $1');
labels=regexprep(labels,'Fornix.*','Fornix / stria terminalis');
labels=strrep(labels,'Superior fronto-occipital fasciculus','SFOF');
labels=strrep(labels,'Retrolenticular part of internal capsule','Retrolenticular IC');
labels=strrep(labels,'Anterior limb of internal capsule','Anterior limb IC');
labels=strrep(labels,'Posterior limb of internal capsule','Posterior limb IC');
labels=strrep(labels,'Posterior thalamic radiation','Posterior thalamic rad.');
labels=strrep(labels,'Superior corona radiata','Sup. corona radiata');
labels=strrep(labels,'Anterior corona radiata','Ant. corona radiata');
labels=strrep(labels,'Cingulum (cingulate gyrus)','Cingulum');
labels=strrep(labels,'External capsule','External capsule');
end
