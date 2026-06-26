function fig_atlas_wm_panel(atlas_nii,wm_nii,out_png,title_text,opts)
arguments
    atlas_nii char
    wm_nii char
    out_png char
    title_text char = 'GM atlas and WM target mask'
    opts.cut char = '-8,-20,18'
    opts.python char = ''
end
root=fig_project_root();
helper=fullfile(root,'code','figure_helpers','atlas_wm_panel.py');
python_exe = local_python_executable(opts.python);
cmd=sprintf('"%s" "%s" --atlas "%s" --wm "%s" --output "%s" --title "%s" --cut "%s"', ...
    python_exe,helper,atlas_nii,wm_nii,out_png,title_text,opts.cut);
try
    cfg=load_project_config();
    bg_file=fullfile(cfg.data_dir,'T1_avg.nii');
    if exist(bg_file,'file'), cmd=sprintf('%s --bg "%s"',cmd,bg_file); end
catch
end
[status,msg]=system(cmd);
if status~=0
    error('Atlas-WM panel failed:\n%s\n%s',cmd,msg);
end
end

function python_exe = local_python_executable(preferred)
candidates = {};
if strlength(string(preferred)) > 0, candidates{end+1} = char(preferred); end
env_py = getenv('RIDGE_GM_WM_PYTHON');
if ~isempty(env_py), candidates{end+1} = env_py; end
candidates{end+1} = fullfile(getenv('USERPROFILE'),'AppData','Local','Programs','Python','Python314','python.exe');
candidates{end+1} = fullfile(getenv('USERPROFILE'),'AppData','Local','Programs','Python','Python313','python.exe');
candidates{end+1} = fullfile(getenv('USERPROFILE'),'AppData','Local','Programs','Python','Python312','python.exe');
candidates{end+1} = 'python';
candidates{end+1} = 'py';
for i = 1:numel(candidates)
    c = candidates{i};
    if any(c == filesep) || contains(c,':')
        if exist(c,'file'), python_exe = c; return; end
    else
        python_exe = c; return;
    end
end
python_exe = 'python';
end
