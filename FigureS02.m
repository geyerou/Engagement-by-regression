%% FigureS02 - Ridge parameters and model summary
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS02');
L=readtable(fullfile(cfg.result_dirs{5},'shared_lambda_refit','shared_lambda_record.csv'),'Delimiter',',','TextType','string');
fig=figure('Position',[100 100 330 250]); histogram(log10(L.subject_final_lambda),18,'FaceColor',[.55 .55 .55],'EdgeColor','none'); hold on; xline(log10(L.shared_lambda(1)),'Color',[.85 .25 .20],'LineWidth',1.2);
xlabel('log10 selected lambda'); ylabel('Subjects'); title('Subject final lambdas'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_a_subject_lambda_distribution.png'));
fig=figure('Position',[100 100 330 250]); boxchart(log10(L.subject_final_lambda),'BoxFaceColor',[.70 .70 .70],'BoxEdgeColor',[.35 .35 .35],'MarkerStyle','.'); hold on; yline(log10(L.shared_lambda(1)),'Color',[.85 .25 .20],'LineWidth',1.2); ylabel('log10 lambda'); xticklabels({'Subjects'}); title('Shared lambda reference'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_b_shared_lambda_reference.png'));
M=readtable(fullfile(cfg.result_dirs{19},'model_comparison_table.csv'),'Delimiter',',','TextType','string');
G=groupsummary(M,'model','mean','mean_R2');
fig=figure('Position',[100 100 360 260]); labels=strrep(string(G.model),'_',' '); bar(categorical(labels),G.mean_mean_R2,'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylabel('Mean subject R²'); xtickangle(20); title('Model mean performance'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_c_model_mean_performance.png'));
P=groupsummary(M,'model','mean','positive_R2_fraction');
fig=figure('Position',[100 100 360 260]); labels=strrep(string(P.model),'_',' '); bar(categorical(labels),P.mean_positive_R2_fraction,'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylabel('Mean positive R² fraction'); xtickangle(20); title('Positive prediction extent'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_d_model_positive_extent.png'));
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S2 layout guide\n\nArrange lambda distribution, shared-lambda reference, mean model performance, and positive prediction extent in a 2x2 grid.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S2 legend\n\nRidge-model parameter diagnostics. (a) Subject-level final lambda distribution. (b) Shared lambda used for beta-profile refit shown relative to subject lambdas. (c-d) Model performance and positive prediction extent across model variants.\n");
fprintf('FigureS02 complete\n');
