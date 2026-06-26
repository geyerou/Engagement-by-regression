%% s005_write_r2_maps
clear; clc;
cfg = load_project_config();
in_dir = cfg.result_dirs{5};
out_dir = cfg.result_dirs{6};
wm = niftiread(cfg.wm_mask_final) > 0;
all_r2 = zeros(numel(cfg.subject_list),nnz(wm),'single');

for s = 1:numel(cfg.subject_list)
    subject_id = cfg.subject_list{s};
    f = fullfile(in_dir,sprintf('sub-%s_ridge_main_results.mat',subject_id));
    D = load(f,'R2_cv_voxel','prediction_correlation','wm_voxel_indices');
    write_map_like(D.R2_cv_voxel,D.wm_voxel_indices,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('sub-%s_ridge_R2_raw_map.nii.gz',subject_id)));
    write_map_like(max(D.R2_cv_voxel,0),D.wm_voxel_indices,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('sub-%s_ridge_R2_positive_map.nii.gz',subject_id)));
    write_map_like(D.prediction_correlation,D.wm_voxel_indices,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('sub-%s_prediction_correlation_map.nii.gz',subject_id)));
    all_r2(s,:) = D.R2_cv_voxel;
end
write_map_like(mean(all_r2,1),D.wm_voxel_indices,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_mean_ridge_R2_raw_map.nii.gz'));
write_map_like(mean(max(all_r2,0),1),D.wm_voxel_indices,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_mean_ridge_R2_positive_map.nii.gz'));
