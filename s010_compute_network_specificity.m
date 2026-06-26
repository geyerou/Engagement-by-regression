%% s010_compute_network_specificity
% Uses positive Delta R2 only. Absolute negative ablations are not evidence
% of contribution and would artificially inflate specificity.
clear; clc;
cfg=load_project_config();
in_dir=cfg.result_dirs{9}; out_dir=cfg.result_dirs{11};

for s=1:numel(cfg.subject_list)
    subject_id=cfg.subject_list{s};
    for scheme="17net"
        D=load(fullfile(in_dir,sprintf( ...
            'sub-%s_network_contribution_%s.mat',subject_id,scheme)), ...
            'Delta_R2','wm_voxel_indices');
        C=max(single(D.Delta_R2),0);
        total=sum(C,2);
        p=C./max(total,eps('single'));
        H=-sum(p.*log(max(p,eps('single'))),2);
        H_norm=H/log(size(C,2));
        Specificity=1-H_norm;
        Specificity(total<=0)=NaN;
        save(fullfile(out_dir,sprintf( ...
            'sub-%s_network_specificity_%s.mat',subject_id,scheme)), ...
            'H','H_norm','Specificity','C','-v7.3');
        write_map_like(H_norm,D.wm_voxel_indices,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf('sub-%s_%s_network_entropy_map.nii.gz', ...
            subject_id,scheme)));
        write_map_like(Specificity,D.wm_voxel_indices,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf( ...
            'sub-%s_%s_network_specificity_map.nii.gz',subject_id,scheme)));
    end
end
