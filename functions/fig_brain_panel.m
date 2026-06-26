function fig_brain_panel(input_nii,out_png,title_text,cmap,opts)
arguments
    input_nii char
    out_png char
    title_text char = ''
    cmap char = 'viridis'
    opts.vmin double = NaN
    opts.vmax double = NaN
    opts.threshold double = NaN
    opts.kind char = 'stat'
    opts.symmetric logical = false
    opts.display char = 'ortho'
    opts.cut char = '-8,-20,18'
    opts.python char = ''
end
root=fig_project_root();
helper=fullfile(root,'code','figure_helpers','brain_panel.py');
python_exe = local_python_executable(opts.python);
cmd=sprintf('"%s" "%s" --input "%s" --output "%s" --title "%s" --cmap "%s"', ...
    python_exe, ...
    helper,input_nii,out_png,title_text,cmap);
try
    cfg=load_project_config();
    bg_file=fullfile(cfg.data_dir,'T1_avg.nii');
    if exist(bg_file,'file')
        cmd=sprintf('%s --bg "%s" --display "%s" --cut "%s"',cmd,bg_file,opts.display,opts.cut);
    end
catch
end
if ~isnan(opts.vmin), cmd=sprintf('%s --vmin %.12g',cmd,opts.vmin); end
if ~isnan(opts.vmax), cmd=sprintf('%s --vmax %.12g',cmd,opts.vmax); end
if ~isnan(opts.threshold), cmd=sprintf('%s --threshold %.12g',cmd,opts.threshold); end
if opts.symmetric, cmd=sprintf('%s --symmetric',cmd); end
cmd=sprintf('%s --kind "%s"',cmd,opts.kind);
[status,msg]=system(cmd);
if status~=0
    error('Brain panel failed:\n%s\n%s',cmd,msg);
end
end

function python_exe = local_python_executable(preferred)
candidates = {};
if strlength(string(preferred)) > 0
    candidates{end+1} = char(preferred);
end
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
        if exist(c,'file')
            python_exe = c;
            return;
        end
    else
        python_exe = c;
        return;
    end
end
python_exe = 'python';
end
