function fig_paired_cloud(vals,labels,ylabel_text,title_text,out_png)
% Compact paired raincloud-like plot for two related measures.
vals = double(vals);
fig=figure('Position',[100 100 330 270]); hold on;
C=fig_colors();
for j=1:2
    x0=j;
    v=vals(:,j); v=v(isfinite(v));
    if numel(v)>2
        [f,xi]=ksdensity(v,'NumPoints',80);
        f=f./max(f)*0.18;
        patch(x0 - f,xi,C.lightgray,'EdgeColor','none','FaceAlpha',0.9);
    end
end
for i=1:size(vals,1)
    plot([1 2],vals(i,:),'-','Color',[.78 .78 .78 .35],'LineWidth',0.55);
end
jit=(rand(size(vals))-0.5)*0.08;
scatter(ones(size(vals,1),1)+jit(:,1),vals(:,1),12,'filled', ...
    'MarkerFaceColor',[.35 .35 .35],'MarkerFaceAlpha',0.45,'MarkerEdgeColor','none');
scatter(2*ones(size(vals,1),1)+jit(:,2),vals(:,2),12,'filled', ...
    'MarkerFaceColor',[.35 .35 .35],'MarkerFaceAlpha',0.45,'MarkerEdgeColor','none');
boxchart([ones(size(vals,1),1);2*ones(size(vals,1),1)],[vals(:,1);vals(:,2)], ...
    'BoxFaceColor','none','BoxEdgeColor',[.15 .15 .15],'WhiskerLineColor',[.15 .15 .15], ...
    'MarkerStyle','none','LineWidth',0.8);
xlim([0.55 2.35]); xticks([1 2]); xticklabels(labels);
ylabel(ylabel_text); title(title_text); set(gca,'TickDir','out');
fig_export(fig,out_png);
end
