%% s006_extract_beta_profiles
clear; clc;
cfg = load_project_config();
labels = read_schaefer_labels(cfg);
in_dir = cfg.result_dirs{5};
out_dir = cfg.result_dirs{7};
shared_dir = fullfile(in_dir,'shared_lambda_refit');

for s = 1:numel(cfg.subject_list)
    subject_id = cfg.subject_list{s};
    switch lower(cfg.beta_source_for_downstream)
        case 'subject_specific'
            D = load(fullfile(in_dir,sprintf( ...
                'sub-%s_ridge_main_results.mat',subject_id)), ...
                'B_ridge_final','wm_voxel_indices','final_lambda');
            Beta_profile = D.B_ridge_final';
            wm_voxel_indices = D.wm_voxel_indices;
            best_lambda = D.final_lambda;
            beta_source = 'subject_specific';
        case 'shared_lambda'
            D = load(fullfile(shared_dir,sprintf( ...
                'sub-%s_shared_lambda_beta.mat',subject_id)), ...
                'B_ridge_shared','wm_voxel_indices','shared_lambda');
            Beta_profile = D.B_ridge_shared';
            wm_voxel_indices = D.wm_voxel_indices;
            best_lambda = D.shared_lambda;
            beta_source = 'shared_lambda';
        otherwise
            error('Unknown beta_source_for_downstream: %s', ...
                cfg.beta_source_for_downstream);
    end
    roi_labels = labels.roi_names;
    network_labels_17 = labels.network17_id;
    save(fullfile(out_dir,sprintf('sub-%s_beta_profile.mat',subject_id)), ...
        'Beta_profile','roi_labels','network_labels_17', ...
        'wm_voxel_indices','best_lambda', ...
        'beta_source','-v7.3');
end
