%% s014_null_model_ridge_prediction
% Exact circular-shift null, accumulated in resumable checkpoints.
% One invocation automatically continues each subject to cfg.n_null.
clear; clc;
cfg=load_project_config();
if ~cfg.do_null_model, fprintf('Null model disabled.\n'); return; end
if ~strcmpi(cfg.null_type,'circular_shift')
    error('Implemented primary null_type is circular_shift.');
end
rng('shuffle');
input_dir=cfg.result_dirs{4}; ridge_dir=cfg.result_dirs{5};
out_dir=cfg.result_dirs{15};
min_shift=round(cfg.min_shift_seconds/cfg.TR);

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    state_file=fullfile(out_dir,sprintf( ...
        'sub-%s_null_circular_shift_state.mat',subject_id));
    observed=load(fullfile(ridge_dir,sprintf( ...
        'sub-%s_ridge_main_results.mat',subject_id)), ...
        'R2_cv_voxel','best_lambda_outer','wm_voxel_indices');
    if exist(state_file,'file')
        state=load(state_file);
        exceed_count=state.exceed_count;
        n_done=state.n_done;
    else
        exceed_count=zeros(size(observed.R2_cv_voxel),'uint16');
        n_done=0;
    end
    if n_done>=cfg.n_null
        fprintf('Null %s already complete (%d).\n',subject_id,n_done);
        continue;
    end

    fprintf('\nNull %s: resuming at %d/%d\n', ...
        subject_id,n_done,cfg.n_null);
    [X_runs,Y_runs,~]=load_subject_runs(fullfile(input_dir, ...
        sprintf('sub-%s_model_input.mat',subject_id)));
    n_runs=numel(X_runs);
    models=cell(n_runs,1);
    train_indices=cell(n_runs,1);
    for test_run=1:n_runs
        train_indices{test_run}=setdiff(1:n_runs,test_run);
        models{test_run}=prepare_ridge_fold( ...
            vertcat(X_runs{train_indices{test_run}}),X_runs{test_run}, ...
            vertcat(Y_runs{train_indices{test_run}}), ...
            observed.best_lambda_outer(test_run));
    end

    while n_done<cfg.n_null
        n_add=min(cfg.null_batch_size,cfg.n_null-n_done);
        fprintf('  batch %d -> %d\n',n_done,n_done+n_add);
        for b=1:n_add
            Ynull=cellfun(@(y)random_circular_shift(y,min_shift), ...
                Y_runs,'UniformOutput',false);
            sse=zeros(size(observed.R2_cv_voxel),'double');
            sst=zeros(size(observed.R2_cv_voxel),'double');
            for test_run=1:n_runs
                pred=predict_prepared_ridge(models{test_run}, ...
                    vertcat(Ynull{train_indices{test_run}}), ...
                    cfg.voxel_block_size);
                y=double(Ynull{test_run});
                sse=sse+sum((y-double(pred)).^2,1);
                sst=sst+sum((y-mean(y,1)).^2,1);
            end
            null_r2=single(1-sse./max(sst,eps));
            exceed_count=exceed_count+uint16( ...
                null_r2>=observed.R2_cv_voxel);
        end
        n_done=n_done+n_add;
        p_voxel=single((1+double(exceed_count))/(1+n_done));
        q_voxel=single(fdr_bh(double(p_voxel)));
        save(state_file,'exceed_count','n_done','p_voxel','q_voxel', ...
            'min_shift','-v7.3');
        write_map_like(p_voxel,observed.wm_voxel_indices, ...
            cfg.wm_mask_final,fullfile(out_dir,sprintf( ...
            'sub-%s_null_pmap_R2.nii.gz',subject_id)));
        write_map_like(q_voxel,observed.wm_voxel_indices, ...
            cfg.wm_mask_final,fullfile(out_dir,sprintf( ...
            'sub-%s_null_qmap_R2.nii.gz',subject_id)));
    end
    clear X_runs Y_runs models
end
