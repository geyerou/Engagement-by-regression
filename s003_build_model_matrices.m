%% s003_build_model_matrices
% Build lightweight manifests. Do not duplicate the large Y matrices.
clear; clc;
cfg = load_project_config();
labels = read_schaefer_labels(cfg);
in_dir = cfg.result_dirs{3};
out_dir = cfg.result_dirs{4};

for s = 1:numel(cfg.subject_list)
    subject_id = cfg.subject_list{s};
    run_files = cell(numel(cfg.run_ids),1);
    n_timepoints = zeros(numel(cfg.run_ids),1);
    wm_voxel_indices = [];
    for r = 1:numel(cfg.run_ids)
        run_files{r} = fullfile(in_dir, sprintf( ...
            'sub-%s_run-%s_timeseries.mat',subject_id,cfg.run_ids{r}));
        if ~exist(run_files{r},'file')
            error('Missing timeseries file: %s',run_files{r});
        end
        m = matfile(run_files{r});
        sz = size(m,'X_gm_raw');
        n_timepoints(r) = sz(1);
        idx = m.wm_voxel_indices;
        if isempty(wm_voxel_indices)
            wm_voxel_indices = idx;
        elseif ~isequal(wm_voxel_indices,idx)
            error('WM indices differ across runs for subject %s.',subject_id);
        end
    end
    model_input = struct();
    model_input.subject_id = subject_id;
    model_input.run_ids = cfg.run_ids;
    model_input.run_files = run_files;
    model_input.n_timepoints = n_timepoints;
    model_input.TR = cfg.TR;
    model_input.wm_voxel_indices = wm_voxel_indices;
    model_input.roi_labels = labels.roi_names;
    model_input.network_labels_17 = labels.network17_id;
    model_input.network_names_17 = labels.network17_names;
    save(fullfile(out_dir,sprintf('sub-%s_model_input.mat',subject_id)), ...
        '-struct','model_input','-v7.3');
end
