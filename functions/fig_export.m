function fig_export(fig,out_png)
fig_apply_publication_style(fig);
set(fig,'Color','w','Units','centimeters');
exportgraphics(fig,out_png,'Resolution',600,'BackgroundColor','white');
[p,n]=fileparts(out_png);
try
    exportgraphics(fig,fullfile(p,[n '.pdf']),'ContentType','vector', ...
        'BackgroundColor','white');
catch
end
savefig(fig,fullfile(p,[n '.fig']));
close(fig);
end
