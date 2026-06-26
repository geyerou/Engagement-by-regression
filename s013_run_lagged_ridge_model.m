%% s013_run_lagged_ridge_model
% Computationally intensive. Test one subject first.
clear; clc;
cfg=load_project_config();
if ~cfg.do_lagged_model, fprintf('Lagged model disabled.\n'); return; end
labels=read_schaefer_labels(cfg);
input_dir=cfg.result_dirs{4};
out_dir=cfg.result_dirs{14}; lags=cfg.lag_list_main;

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    fprintf('\nLagged model %s (%d/%d)\n',subject_id,s,numel(cfg.subject_list));
    [X_runs,Y_runs,meta]=load_subject_runs(fullfile(input_dir, ...
        sprintf('sub-%s_model_input.mat',subject_id)));
    if strcmpi(cfg.lagged_feature_space,'network17')
        X_runs=cellfun(@(x)aggregate_network_timeseries(x, ...
            labels.network17_id),X_runs,'UniformOutput',false);
    end
    n_features=size(X_runs{1},2);
    XL=cell(size(X_runs)); X0=cell(size(X_runs)); YL=cell(size(Y_runs));
    max_abs=max(abs(lags));
    for r=1:numel(X_runs)
        [XL{r},YL{r}]=build_lagged_run(X_runs{r},Y_runs{r},lags);
        valid=(1+max_abs):(size(X_runs{r},1)-max_abs);
        X0{r}=X_runs{r}(valid,:);
    end
    result=crossval_nested_ridge_r2(XL,YL,cfg.lagged_lambda_grid, ...
        cfg.lagged_lambda_selection_voxels,cfg.voxel_block_size);
    % Refit zero-lag on the identical trimmed time points. Comparing against
    % the full-length main model would mix model effects with sample changes.
    zero_lag_matched=crossval_nested_ridge_r2(X0,YL,cfg.lambda_grid, ...
        cfg.lagged_lambda_selection_voxels,cfg.voxel_block_size);
    Delta_R2_lagged=result.R2-zero_lag_matched.R2;
    final_lambda=exp(mean(log(result.best_lambda_outer)));
    [lag_energy,best_lag]=fit_lag_beta_summary(vertcat(XL{:}), ...
        vertcat(YL{:}),final_lambda,cfg.voxel_block_size,n_features, ...
        numel(lags),lags);
    wm_voxel_indices=meta.wm_voxel_indices;
    feature_space=cfg.lagged_feature_space;
    save(fullfile(out_dir,sprintf('sub-%s_lagged_ridge_results.mat',subject_id)), ...
        'result','zero_lag_matched','Delta_R2_lagged','lag_energy','best_lag','lags', ...
        'feature_space','wm_voxel_indices','-v7.3');
    write_map_like(result.R2,wm_voxel_indices,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('sub-%s_lagged_ridge_R2_map.nii.gz',subject_id)));
    write_map_like(Delta_R2_lagged,wm_voxel_indices,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf( ...
        'sub-%s_lagged_minus_zerolag_R2_map.nii.gz',subject_id)));
    write_map_like(best_lag,wm_voxel_indices,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('sub-%s_best_lag_map.nii.gz',subject_id)));
    clear X_runs Y_runs XL X0 YL
end
