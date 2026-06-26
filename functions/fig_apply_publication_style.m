function fig_apply_publication_style(fig)
% Apply a compact journal-style finish before export.
if nargin < 1 || isempty(fig), fig = gcf; end
set(fig,'Color','w','InvertHardcopy','off');
axs = findall(fig,'Type','axes');
for i = 1:numel(axs)
    ax = axs(i);
    tag = string(get(ax,'Tag'));
    if tag == "legend"
        continue;
    end
    set(ax,'FontName','Arial','FontSize',10,'LineWidth',0.55, ...
        'TickDir','out','Layer','top','Box','off');
    ax.Title.String = '';
    ax.Title.FontWeight = 'normal';
    ax.Title.FontSize = 10;
    ax.XLabel.FontSize = 10;
    ax.YLabel.FontSize = 10;
    try
        ax.Toolbar.Visible = 'off';
    catch
    end
end
legs = findall(fig,'Type','Legend');
for i = 1:numel(legs)
    set(legs(i),'Box','off','FontName','Arial','FontSize',9);
end
cbs = findall(fig,'Type','ColorBar');
for i = 1:numel(cbs)
    set(cbs(i),'FontName','Arial','FontSize',9,'LineWidth',0.45);
    try
        cbs(i).Ticks = linspace(cbs(i).Limits(1),cbs(i).Limits(2),4);
    catch
    end
end
bars = findall(fig,'Type','Bar');
for i = 1:numel(bars)
    try
        bars(i).BarWidth = 0.52;
    catch
    end
end
texts = findall(fig,'Type','Text');
for i = 1:numel(texts)
    set(texts(i),'FontName','Arial');
end
end
