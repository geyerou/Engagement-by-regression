%% s011_compute_fc_baseline
clear; clc;
cfg=load_project_config();
input_dir=cfg.result_dirs{4}; out_dir=cfg.result_dirs{12};

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    [X_runs,Y_runs,meta]=load_subject_runs(fullfile(input_dir, ...
        sprintf('sub-%s_model_input.mat',subject_id)));
    X=vertcat(X_runs{:});
    X=X-mean(X,1); X=X./max(std(X,0,1),eps('single'));
    n_vox=size(Y_runs{1},2);
    FC_profile=zeros(n_vox,size(X,2),'single');
    Y=vertcat(Y_runs{:});
    for first=1:cfg.voxel_block_size:n_vox
        idx=first:min(first+cfg.voxel_block_size-1,n_vox);
        y=Y(:,idx); y=y-mean(y,1); y=y./max(std(y,0,1),eps('single'));
        FC_profile(idx,:)=single((double(y)'*double(X))/(size(X,1)-1));
    end
    wm_voxel_indices=meta.wm_voxel_indices;
    save(fullfile(out_dir,sprintf('sub-%s_FC_profile.mat',subject_id)), ...
        'FC_profile','wm_voxel_indices','-v7.3');
    clear X Y X_runs Y_runs FC_profile
end
