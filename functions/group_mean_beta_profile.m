function [group_profile,wm_idx] = group_mean_beta_profile(cfg)
in_dir=cfg.result_dirs{7};
acc=[]; wm_idx=[];
for s=1:numel(cfg.subject_list)
    D=load(fullfile(in_dir,sprintf('sub-%s_beta_profile.mat', ...
        cfg.subject_list{s})),'Beta_profile','wm_voxel_indices');
    if isempty(acc)
        acc=zeros(size(D.Beta_profile),'double');
        wm_idx=D.wm_voxel_indices;
    end
    acc=acc+double(D.Beta_profile);
end
group_profile=single(acc/numel(cfg.subject_list));
end
