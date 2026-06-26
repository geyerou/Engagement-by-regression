%% s002_extract_timeseries_with_wm_smoothing
% Processes subjects serially. Each 4D run is cleared before the next run.
clear; clc;
cfg = load_project_config();
start_optional_thread_pool(cfg);
out_dir = cfg.result_dirs{3};

wm_file = cfg.wm_mask_final;
if ~exist(wm_file,'file')
    error('Final WM mask missing. Run s001 first.');
end

for s = 1:numel(cfg.subject_list)
    subject_id = cfg.subject_list{s};
    fprintf('\nSubject %s (%d/%d)\n',subject_id,s,numel(cfg.subject_list));
    for r = 1:numel(cfg.run_dirs)
        out_file = fullfile(out_dir, sprintf('sub-%s_run-%s_timeseries.mat', ...
            subject_id, cfg.run_ids{r}));
        if exist(out_file,'file')
            fprintf('  %s exists; skip.\n',cfg.run_ids{r});
            continue;
        end
        func_file = get_run_file(cfg,subject_id,r);
        fprintf('  Extracting %s\n',cfg.run_ids{r});
        run_data = extract_run_timeseries(func_file, cfg.schaefer400_file, ...
            wm_file, cfg);
        run_data.subject_id = subject_id;
        run_data.run_id = cfg.run_ids{r};
        save(out_file, '-struct', 'run_data', '-v7.3');
        clear run_data
    end
end
