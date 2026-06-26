%% s008_compute_network_contribution
% Leave-one-network-out predictive ablation using observed outer-fold lambdas.
clear; clc;
cfg=load_project_config();
labels=read_schaefer_labels(cfg);
input_dir=cfg.result_dirs{4}; ridge_dir=cfg.result_dirs{5};
out_dir=cfg.result_dirs{9};

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    fprintf('\nNetwork ablation %s (%d/%d)\n',subject_id,s,numel(cfg.subject_list));
    [X_runs,Y_runs,meta]=load_subject_runs(fullfile(input_dir, ...
        sprintf('sub-%s_model_input.mat',subject_id)));
    observed=load(fullfile(ridge_dir,sprintf( ...
        'sub-%s_ridge_main_results.mat',subject_id)), ...
        'R2_cv_voxel','best_lambda_outer');

    net_id=labels.network17_id;
    scheme_name='17net';
    n_net=max(net_id);
    Delta_R2=zeros(size(Y_runs{1},2),n_net,'single');
    R2_without=zeros(size(Y_runs{1},2),n_net,'single');
    for k=1:n_net
        fprintf('  %s remove network %d/%d\n',scheme_name,k,n_net);
        keep=net_id~=k;
        X_reduced=cellfun(@(x)x(:,keep),X_runs,'UniformOutput',false);
        reduced=crossval_fixed_lambda_r2(X_reduced,Y_runs, ...
            observed.best_lambda_outer,cfg.voxel_block_size);
        R2_without(:,k)=reduced(:);
        Delta_R2(:,k)=observed.R2_cv_voxel(:)-reduced(:);
    end
    wm_voxel_indices=meta.wm_voxel_indices;
    network_names=labels.network17_names;
    save(fullfile(out_dir,sprintf( ...
        'sub-%s_network_contribution_%s.mat',subject_id,scheme_name)), ...
        'Delta_R2','R2_without','wm_voxel_indices','network_names','-v7.3');
    if cfg.write_subject_network_maps
        for k=1:n_net
            write_map_like(Delta_R2(:,k),wm_voxel_indices,cfg.wm_mask_final, ...
                fullfile(out_dir,sprintf( ...
                'sub-%s_%s-network-%02d_DeltaR2_map.nii.gz', ...
                subject_id,scheme_name,k)));
        end
    end
    clear X_runs Y_runs
end
