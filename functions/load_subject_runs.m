function [X_runs,Y_runs,meta] = load_subject_runs(model_input_file)
meta = load(model_input_file);
n = numel(meta.run_files);
X_runs = cell(n,1);
Y_runs = cell(n,1);
for r = 1:n
    tmp = load(meta.run_files{r}, 'X_gm_raw','Y_wm_raw');
    X_runs{r} = single(tmp.X_gm_raw);
    Y_runs{r} = single(tmp.Y_wm_raw);
end
end
