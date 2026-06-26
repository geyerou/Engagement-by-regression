%% s016_group_level_analysis
clear; clc;
cfg=load_project_config();
ridge_dir=cfg.result_dirs{5}; net_dir=cfg.result_dirs{9};
spec_dir=cfg.result_dirs{11}; out_dir=cfg.result_dirs{17};
n_sub=numel(cfg.subject_list);
wm=niftiread(cfg.wm_mask_final)>0;
R=zeros(n_sub,nnz(wm),'single'); wm_idx=find(wm);

for s=1:n_sub
    subject_id=cfg.subject_list{s};
    D=load(fullfile(ridge_dir,sprintf( ...
        'sub-%s_ridge_main_results.mat',subject_id)), ...
        'R2_cv_voxel','wm_voxel_indices');
    if ~isequal(wm_idx,D.wm_voxel_indices)
        error('WM indices mismatch for subject %s.',subject_id);
    end
    R(s,:)=D.R2_cv_voxel;
end
group_mean=mean(R,1);
group_median=median(R,1);
group_sd=std(R,0,1);
group_t=group_mean./max(group_sd/sqrt(n_sub),eps('single'));
if license('test','statistics_toolbox')
    p=2*tcdf(-abs(double(group_t)),n_sub-1);
else
    warning('Statistics Toolbox unavailable; using normal approximation.');
    p=erfc(abs(double(group_t))/sqrt(2));
end
q=fdr_bh(p);
consistency=mean(R>0,1);

write_map_like(group_mean,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_mean_R2_map.nii.gz'));
write_map_like(group_median,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_median_R2_map.nii.gz'));
write_map_like(group_t,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_t_R2_map.nii.gz'));
write_map_like(single(q),wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_q_R2_map.nii.gz'));
write_map_like(single(q<0.05 & group_mean>0),wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'significant_R2_mask_FDR.nii.gz'));
write_map_like(consistency,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_positive_R2_consistency_map.nii.gz'));
save(fullfile(out_dir,'group_R2_statistics.mat'),'R','group_mean', ...
    'group_median','group_t','p','q','consistency','wm_idx','-v7.3');

% Group network contribution and specificity maps.
for scheme="17net"
    n_net=17;
    C=zeros(numel(wm_idx),n_net,n_sub,'single');
    for s=1:n_sub
        subject_id=cfg.subject_list{s};
        D=load(fullfile(net_dir,sprintf( ...
            'sub-%s_network_contribution_%s.mat',subject_id,scheme)), ...
            'Delta_R2');
        C(:,:,s)=D.Delta_R2;
    end
    meanC=mean(C,3);
    for k=1:size(meanC,2)
        write_map_like(meanC(:,k),wm_idx,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf( ...
            'group_mean_%s_network-%02d_DeltaR2.nii.gz',scheme,k)));
    end
    save(fullfile(out_dir,sprintf('group_%s_contribution.mat',scheme)), ...
        'meanC','-v7.3');
end

for scheme="17net"
    S=zeros(n_sub,numel(wm_idx),'single');
    for s=1:n_sub
        D=load(fullfile(spec_dir,sprintf( ...
            'sub-%s_network_specificity_%s.mat',cfg.subject_list{s},scheme)), ...
            'Specificity');
        S(s,:)=D.Specificity;
    end
    write_map_like(mean(S,1,'omitnan'),wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('group_mean_%s_specificity.nii.gz',scheme)));
end
