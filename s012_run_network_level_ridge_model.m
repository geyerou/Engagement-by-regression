%% s012_run_network_level_ridge_model
clear; clc;
cfg=load_project_config();
labels=read_schaefer_labels(cfg);
input_dir=cfg.result_dirs{4}; out_dir=cfg.result_dirs{13};

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    [X_runs,Y_runs,meta]=load_subject_runs(fullfile(input_dir, ...
        sprintf('sub-%s_model_input.mat',subject_id)));
    for scheme=17
        ids=labels.network17_id;
        Xnet=cellfun(@(x)aggregate_network_timeseries(x,ids), ...
            X_runs,'UniformOutput',false);
        result=crossval_nested_ridge_r2(Xnet,Y_runs,cfg.lambda_grid, ...
            cfg.lambda_selection_voxels,cfg.voxel_block_size);
        wm_voxel_indices=meta.wm_voxel_indices;
        save(fullfile(out_dir,sprintf( ...
            'sub-%s_network%d_ridge_results.mat',subject_id,scheme)), ...
            'result','wm_voxel_indices','-v7.3');
        write_map_like(result.R2,wm_voxel_indices,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf( ...
            'sub-%s_network%d_ridge_R2_map.nii.gz',subject_id,scheme)));
    end
    clear X_runs Y_runs
end
