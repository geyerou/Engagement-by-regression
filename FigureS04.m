%% FigureS04 - Lag and reduced-network sensitivity
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS04');
M=readtable(fullfile(cfg.result_dirs{19},'model_comparison_table.csv'),'Delimiter',',','TextType','string');
fig=figure('Position',[100 100 430 280]); hold on; mods=unique(M.model,'stable'); X=[];Y=[]; for i=1:numel(mods), X=[X;i*ones(sum(M.model==mods(i)),1)]; Y=[Y;M.mean_R2(M.model==mods(i))]; end
boxchart(X,Y,'BoxFaceColor',[.70 .70 .70],'BoxEdgeColor',[.35 .35 .35],'MarkerStyle','.'); xticks(1:numel(mods)); xticklabels(strrep(mods,'_',' ')); xtickangle(20); ylabel('Subject mean R²'); title('Model sensitivity'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_a_model_comparison_all.png'));
D=unstack(M(:,{'subject_id','model','mean_R2'}),'mean_R2','model'); fig=figure('Position',[100 100 330 260]); histogram(D.ridge_lagged-D.ridge_400_zero_lag,18,'FaceColor',[.55 .55 .55],'EdgeColor','none'); xline(0,'k--'); xlabel('Lagged - zero-lag mean R²'); ylabel('Subjects'); title('Lag model sensitivity'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_b_lag_delta_distribution.png'));
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S4 layout guide\n\nArrange model comparison and lag-minus-zero-lag distribution side by side.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S4 legend\n\nSensitivity analyses for reduced network-level predictors and lagged models. These are not primary because the main story uses strict WM95 zero-lag prediction.\n");
fprintf('FigureS04 complete\n');
