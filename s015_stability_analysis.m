%% s015_stability_analysis
% Split-half standardized beta-profile reliability: runs 1+2 vs runs 3+4.
clear; clc;
cfg=load_project_config();
input_dir=cfg.result_dirs{4}; ridge_dir=cfg.result_dirs{5};
out_dir=cfg.result_dirs{16};

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    [X_runs,Y_runs,meta]=load_subject_runs(fullfile(input_dir, ...
        sprintf('sub-%s_model_input.mat',subject_id)));
    main=load(fullfile(ridge_dir,sprintf( ...
        'sub-%s_ridge_main_results.mat',subject_id)),'final_lambda');
    [B1,~,~,~,~]=fit_final_beta_blocked(vertcat(X_runs{1:2}), ...
        vertcat(Y_runs{1:2}),main.final_lambda,cfg.voxel_block_size);
    [B2,~,~,~,~]=fit_final_beta_blocked(vertcat(X_runs{3:4}), ...
        vertcat(Y_runs{3:4}),main.final_lambda,cfg.voxel_block_size);
    beta_stability=row_correlation(B1',B2');
    wm_voxel_indices=meta.wm_voxel_indices;
    save(fullfile(out_dir,sprintf('sub-%s_beta_stability.mat',subject_id)), ...
        'beta_stability','wm_voxel_indices','-v7.3');
    write_map_like(beta_stability,wm_voxel_indices,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('sub-%s_beta_stability_map.nii.gz',subject_id)));
    clear X_runs Y_runs B1 B2
end
