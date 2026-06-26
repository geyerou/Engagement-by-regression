%% s007_write_roi_beta_maps
% Group maps are the default; 400 maps x 81 subjects is intentionally avoided.
clear; clc;
cfg=load_project_config();
in_dir=cfg.result_dirs{7}; out_dir=cfg.result_dirs{8};
group_sum=[]; wm_idx=[];

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    D=load(fullfile(in_dir,sprintf('sub-%s_beta_profile.mat',subject_id)), ...
        'Beta_profile','wm_voxel_indices');
    B=D.Beta_profile';
    if isempty(group_sum)
        group_sum=zeros(size(B),'double');
        wm_idx=D.wm_voxel_indices;
    end
    group_sum=group_sum+double(B);
    if cfg.write_subject_roi_beta_maps
        sub_dir=fullfile(out_dir,['sub-' subject_id]);
        ensure_dir(sub_dir);
        for roi=1:400
            write_map_like(B(roi,:),wm_idx,cfg.wm_mask_final, ...
                fullfile(sub_dir,sprintf('sub-%s_ROI-%03d_beta_map.nii.gz', ...
                subject_id,roi)));
        end
    end
end

group_beta=single(group_sum/numel(cfg.subject_list));
save(fullfile(out_dir,'group_mean_roi_beta.mat'),'group_beta','wm_idx','-v7.3');
if cfg.write_group_roi_beta_maps
    map_dir=fullfile(out_dir,'group_ROI_beta_maps'); ensure_dir(map_dir);
    for roi=1:400
        write_map_like(group_beta(roi,:),wm_idx,cfg.wm_mask_final, ...
            fullfile(map_dir,sprintf('group_ROI-%03d_beta_map.nii.gz',roi)));
    end
end
