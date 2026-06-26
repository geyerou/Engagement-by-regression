function fig_density_scatter(x,y,xlab,ylab,stat_text,out_png)
x = double(x(:));
y = double(y(:));
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
fig = figure('Position',[100 100 330 300]);
hold on;
nbins = 80;
histogram2(x,y,[nbins nbins],'DisplayStyle','tile','ShowEmptyBins','off', ...
    'EdgeColor','none');
view(2);
colormap(flipud(gray(256)));
cb = colorbar;
cb.Label.String = 'Voxels';
xlabel(xlab);
ylabel(ylab);
text(0.04,0.94,stat_text,'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','top','FontSize',7,'FontWeight','bold','Color',[0.12 0.12 0.12]);
set(gca,'TickDir','out');
fig_export(fig,out_png);
end
