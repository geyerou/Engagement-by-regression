%% s004_run_ridge_main_model
% Nested leave-one-run-out CV. Subjects are intentionally serial.
clear; clc;
cfg = load_project_config();
start_optional_thread_pool(cfg);
in_dir = cfg.result_dirs{4};
out_dir = cfg.result_dirs{5};
rng(20260620,'twister');

for s = 1:numel(cfg.subject_list)
    subject_id = cfg.subject_list{s};
    out_file = fullfile(out_dir,sprintf('sub-%s_ridge_main_results.mat',subject_id));
    if exist(out_file,'file')
        fprintf('Subject %s exists; skip.\n',subject_id);
        continue;
    end
    fprintf('\nRidge subject %s (%d/%d)\n',subject_id,s,numel(cfg.subject_list));
    input_file = fullfile(in_dir,sprintf('sub-%s_model_input.mat',subject_id));
    [X_runs,Y_runs,meta] = load_subject_runs(input_file);
    n_runs = numel(X_runs);
    n_vox = size(Y_runs{1},2);
    sample_n = min(cfg.lambda_selection_voxels,n_vox);
    voxel_sample = round(linspace(1,n_vox,sample_n));

    sse = zeros(1,n_vox,'double');
    sst = zeros(1,n_vox,'double');
    sum_y = zeros(1,n_vox,'double');
    sum_p = zeros(1,n_vox,'double');
    sum_yy = zeros(1,n_vox,'double');
    sum_pp = zeros(1,n_vox,'double');
    sum_yp = zeros(1,n_vox,'double');
    n_total = 0;
    best_lambda_outer = zeros(n_runs,1);
    lambda_performance = zeros(n_runs,numel(cfg.lambda_grid));
    if cfg.save_cv_predictions
        Y_pred_cv = cell(n_runs,1);
    else
        Y_pred_cv = {};
    end

    for test_run = 1:n_runs
        train_runs = setdiff(1:n_runs,test_run);
        [best_lambda_outer(test_run),lambda_performance(test_run,:)] = ...
            select_lambda_inner_cv(X_runs,Y_runs,train_runs, ...
            cfg.lambda_grid,voxel_sample);
        fprintf('  Outer run %d lambda = %.4g\n', ...
            test_run,best_lambda_outer(test_run));
        Xtr = vertcat(X_runs{train_runs});
        Ytr = vertcat(Y_runs{train_runs});
        Xte = X_runs{test_run};
        Yte = Y_runs{test_run};
        pred = ridge_fit_predict_blocked(Xtr,Ytr,Xte, ...
            best_lambda_outer(test_run),cfg.voxel_block_size,false);
        if cfg.save_cv_predictions, Y_pred_cv{test_run}=pred; end
        yd = double(Yte); pd = double(pred);
        sse = sse + sum((yd-pd).^2,1);
        centered = yd-mean(yd,1);
        sst = sst + sum(centered.^2,1);
        sum_y = sum_y + sum(yd,1); sum_p = sum_p + sum(pd,1);
        sum_yy = sum_yy + sum(yd.^2,1); sum_pp = sum_pp + sum(pd.^2,1);
        sum_yp = sum_yp + sum(yd.*pd,1);
        n_total = n_total + size(yd,1);
        clear Xtr Ytr Xte Yte pred yd pd
    end

    R2_cv_voxel = single(1-sse./max(sst,eps));
    numerator = n_total*sum_yp-sum_y.*sum_p;
    denominator = sqrt(max(n_total*sum_yy-sum_y.^2,0).* ...
                       max(n_total*sum_pp-sum_p.^2,0));
    prediction_correlation = single(numerator./max(denominator,eps));

    final_lambda = exp(mean(log(best_lambda_outer)));
    Xall = vertcat(X_runs{:});
    Yall = vertcat(Y_runs{:});
    [B_ridge_final,muX_final,sdX_final,muY_final,sdY_final] = ...
        fit_final_beta_blocked(Xall,Yall,final_lambda,cfg.voxel_block_size);

    wm_voxel_indices = meta.wm_voxel_indices;
    lambda_grid = cfg.lambda_grid;
    save(out_file,'subject_id','R2_cv_voxel','prediction_correlation', ...
        'B_ridge_final','best_lambda_outer','final_lambda','lambda_grid', ...
        'lambda_performance','wm_voxel_indices','muX_final','sdX_final', ...
        'muY_final','sdY_final','Y_pred_cv','-v7.3');
    clear X_runs Y_runs Xall Yall B_ridge_final
end
