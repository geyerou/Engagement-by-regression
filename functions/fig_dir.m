function [figdir,paneldir]=fig_dir(fig_name)
root=fig_project_root();
figdir=fullfile(root,'figures',fig_name);
paneldir=fullfile(figdir,'panels');
ensure_dir(figdir); ensure_dir(paneldir);
end
