%% FigureS05 - Negative / cautionary analyses
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS05');
K=readtable(fullfile(cfg.output_dir,'result_s009_beta_similarity_gradient','s009b_beta_kmeans_parcellation','beta_kmeans_summary.csv'),'Delimiter',',','TextType','string');
[~,ord]=sort(K.k); K=K(ord,:);
fig=figure('Position',[100 100 320 250]); cats=categorical(string(K.k)); cats=reordercats(cats,string(K.k)); bar(cats,K.mean_sample_silhouette,'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylabel('Mean silhouette'); xlabel('k'); title('Weak k-means parcellation support'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_a_kmeans_silhouette.png'));
B=readtable(fullfile(cfg.result_dirs{18},'behavior_prediction_summary.csv'),'Delimiter',',','TextType','string');
fig=figure('Position',[100 100 430 260]); hold on;
labs=local_outcome_labels(B.outcome); y=1:height(B);
plot([0 0],[0.5 height(B)+0.5],'--','Color',[.55 .55 .55],'LineWidth',.7);
scatter(B.delta_R2,y,34,B.permutation_q_delta_R2,'filled','MarkerEdgeColor',[.20 .20 .20],'LineWidth',.25);
yticks(y); yticklabels(labs); xlabel('Incremental ΔR²'); title('Global behavior prediction'); cb=colorbar; cb.Label.String='Corrected q'; colormap(flipud(gray(256))); set(gca,'TickDir','out','YDir','reverse'); fig_export(fig,fullfile(paneldir,'panel_b_global_behavior_prediction.png'));
LB=readtable(fullfile(cfg.output_dir,'result_s026_behavior_local_prediction','behavior_local_prediction_summary.csv'),'Delimiter',',','TextType','string');
fig=figure('Position',[100 100 500 330]); hold on;
labs=local_local_labels(LB.feature_set,LB.outcome); y=1:height(LB);
plot([0 0],[0.5 height(LB)+0.5],'--','Color',[.55 .55 .55],'LineWidth',.7);
g=double(categorical(LB.feature_set));
colors=[.72 .20 .17; 0 .48 .42; .12 .31 .54];
for k=1:max(g)
    idx=g==k;
    scatter(LB.delta_R2(idx),y(idx),30,'filled','MarkerFaceColor',colors(k,:), ...
        'MarkerEdgeColor',[.20 .20 .20],'LineWidth',.25);
end
yticks(y); yticklabels(labs); xlabel('Incremental ΔR²'); title('Local behavior prediction'); set(gca,'TickDir','out','YDir','reverse');
fig_export(fig,fullfile(paneldir,'panel_c_local_behavior_prediction.png'));
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S5 layout guide\n\nArrange k-means, global behavior, and local behavior negative/cautionary analyses in one row.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S5 legend\n\nAnalyses treated as cautionary or negative: k-means parcellation has weak silhouette support, and global/local behavior prediction did not provide robust corrected associations.\n");
fprintf('FigureS05 complete\n');

function labs=local_outcome_labels(outcome)
labs=strrep(string(outcome),'CogTotalComp_AgeAdj','Total cognition');
labs=strrep(labs,'CogFluidComp_AgeAdj','Fluid cognition');
labs=strrep(labs,'CogCrystalComp_AgeAdj','Crystal cognition');
labs=strrep(labs,'PMAT24_A_CR','PMAT');
end

function labs=local_local_labels(feature_set,outcome)
o=local_outcome_labels(outcome);
f=strrep(string(feature_set),'R2_residual','R2 residual');
f=strrep(f,'beta_minus_FC_G1','Beta-FC G1');
f=strrep(f,'combined','Combined');
labs=f+" / "+o;
end
