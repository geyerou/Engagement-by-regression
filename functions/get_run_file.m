function path_value = get_run_file(cfg, subject_id, run_index)
path_value = fullfile(cfg.data_dir, cfg.run_dirs{run_index}, subject_id, ...
    [cfg.run_ids{run_index} cfg.func_suffix]);
end
