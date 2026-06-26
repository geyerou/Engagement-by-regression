%% s000_config_project
% Create the central configuration. Edit only this file for normal use.

clear; clc;

project_dir = fileparts(fileparts(mfilename('fullpath')));
code_dir = fullfile(project_dir, 'code');
functions_dir = fullfile(code_dir, 'functions');
addpath(functions_dir);

cfg = struct();
cfg.project_dir = project_dir;
cfg.code_dir = code_dir;
cfg.functions_dir = functions_dir;
cfg.data_root = 'F:\Demo_Data_HCPY';
cfg.data_dir = fullfile(cfg.data_root, 'HCP_young_100subs_2ses_2run');
% The project has one active analysis: 4 mm mask-normalized WM smoothing.
cfg.analysis_label = 'wmSmooth4mm';
cfg.output_dir = fullfile(project_dir, 'derivatives');
cfg.docs_dir = fullfile(project_dir, 'docs');

cfg.run_dirs = {'Ses1_LR','Ses1_RL','Ses2_LR','Ses2_RL'};
cfg.run_ids = {'rfMRI_REST1_LR','rfMRI_REST1_RL', ...
               'rfMRI_REST2_LR','rfMRI_REST2_RL'};
cfg.TR = 0.72;
cfg.expected_timepoints = 1200;
cfg.qc_header_subjects_per_run = 1;

cfg.func_suffix = '_preprocessed.nii.gz';
cfg.schaefer400_file = fullfile(cfg.data_dir, 'schaefer_2018', ...
    'Schaefer2018_400_masked_GM05.nii.gz');
cfg.schaefer_order_file = fullfile(cfg.data_dir, 'schaefer_2018', ...
    'Schaefer2018_400Parcels_17Networks_order.txt');
cfg.schaefer_label_file = fullfile(cfg.data_dir, 'schaefer_2018', ...
    'Schaefer2018_400Parcels_17Networks_labels.csv');
cfg.network_schemes = 17;
cfg.wm_mask_source = fullfile(cfg.data_dir, 'Group_WM_Mask_95.nii.gz');
cfg.gm_mask_source = fullfile(cfg.data_dir, 'Group_GM_Mask_thr05.nii.gz');
cfg.motion_table = fullfile(cfg.data_dir, 'HCP_YA_subjects_motion.csv');
cfg.behavior_table = fullfile(cfg.data_dir, 'HCP_YA_subjects.csv');

cfg.result_dirs = cell(19,1);
for k = 0:18
    cfg.result_dirs{k+1} = fullfile(cfg.output_dir, ...
        sprintf('result_s%03d_%s', k, local_stage_name(k)));
    ensure_dir(cfg.result_dirs{k+1});
end
cfg.config_file = fullfile(cfg.result_dirs{1}, 'config.mat');
cfg.wm_mask_final = fullfile(cfg.result_dirs{2}, 'wm_mask_final.nii.gz');

% Main analysis settings.
cfg.main_model = 'ridge';
cfg.cv_scheme = 'nested_leave_one_run_out';
cfg.lambda_grid = logspace(-4, 4, 33);
cfg.lambda_definition = 'XTX_plus_n_times_lambda_I';
cfg.lambda_selection_voxels = 4096;
cfg.voxel_block_size = 2048;
cfg.numeric_class = 'single';
cfg.save_cv_predictions = false; % approximately 0.6 GB/subject if true
cfg.save_fold_betas = false;

% Cross-subject beta comparison.
% Run s004b only after s004 has finished for all configured subjects.
cfg.shared_lambda_method = 'geometric_median';
cfg.beta_source_for_downstream = 'shared_lambda';
% Options: 'subject_specific' or 'shared_lambda'

% Main analysis uses 4 mm FWHM mask-normalized smoothing within WM95.
cfg.do_wm_mask_smoothing = true;
cfg.wm_smooth_fwhm_mm = 4;
cfg.erode_wm_mask = false;
cfg.wm_erosion_radius_vox = 1;

% Parallel policy: never parallelize subjects.
cfg.use_thread_pool = false;
cfg.thread_pool_workers = 4;
cfg.parallelize_voxel_blocks = false;

% Extensions.
cfg.do_lagged_model = true;
cfg.lag_list_main = -5:5;
cfg.lag_list_sensitivity = -7:7;
cfg.lagged_feature_space = 'roi400'; % change to 'network17' for a fast pilot
cfg.lagged_lambda_grid = logspace(-3,3,17);
cfg.lagged_lambda_selection_voxels = 1024;
cfg.do_null_model = true;
cfg.n_null = 1000;
cfg.null_batch_size = 25;
cfg.min_shift_seconds = 20;
cfg.null_type = 'circular_shift';
cfg.do_gradient = true;
cfg.brainspace_path = ...
    'C:\Users\geyerou\Documents\Toolboxes\BrainSpace\matlab';
cfg.brainspace_kernel = 'normalized_angle';
cfg.brainspace_approach = 'dm';
cfg.brainspace_alignment = 'none';
cfg.brainspace_sparsity = 0; % Preserve signed ridge beta fingerprints.
cfg.brainspace_alpha = 0.5;
cfg.brainspace_diffusion_time = 0;
cfg.gradient_components = 10;
cfg.gradient_random_state = 20260620;

cfg.do_beta_parcellation = true;
cfg.parcellation_k = [7 17]; % Number of WM clusters; not Yeo networks.
cfg.parcellation_pca_variance_percent = 95;
cfg.parcellation_pca_max_components = 100;
cfg.parcellation_replicates = 50;
cfg.parcellation_max_iter = 1000;
cfg.parcellation_silhouette_sample = 3000;
cfg.parcellation_random_state = 20260620;
cfg.write_subject_roi_beta_maps = false;
cfg.write_group_roi_beta_maps = true;
cfg.write_subject_network_maps = true;
cfg.fc_beta_top_roi_fraction = 0.05;
cfg.write_subject_fc_beta_maps = true;

% Stability analysis.
cfg.stability_split_names = {'session','balanced'};
% session:  REST1 LR+RL versus REST2 LR+RL
% balanced: REST1 LR+REST2 RL versus REST1 RL+REST2 LR
cfg.stability_use_shared_lambda = true;
cfg.write_subject_stability_maps = true;

% Behavior/cognition analysis.
cfg.behavior_outcomes = { ...
    'CogTotalComp_AgeAdj', ...
    'CogFluidComp_AgeAdj', ...
    'CogCrystalComp_AgeAdj', ...
    'PMAT24_A_CR'};
cfg.behavior_outer_folds = 5;
cfg.behavior_outer_repeats = 20;
cfg.behavior_inner_folds = 5;
cfg.behavior_lambda_grid = logspace(-4,4,25);
cfg.behavior_n_permutations = 1000;
cfg.behavior_random_state = 20260620;

% Safe default: one subject for debugging. Change to Inf only after validation.
cfg.subject_limit = Inf;

subject_dirs = dir(fullfile(cfg.data_dir, cfg.run_dirs{1}));
subject_dirs = subject_dirs([subject_dirs.isdir]);
subject_names = string({subject_dirs.name});
subject_names = subject_names(~startsWith(subject_names,'.'));
subject_names = sort(subject_names);
cfg.all_subject_list = cellstr(subject_names);
cfg.subject_list = cfg.all_subject_list;
if isfinite(cfg.subject_limit)
    cfg.subject_list = cfg.subject_list(1:min(cfg.subject_limit,numel(cfg.subject_list)));
end

ensure_dir(cfg.output_dir);
ensure_dir(cfg.docs_dir);
save(cfg.config_file, 'cfg', '-v7.3');
% Pointer used by all downstream scripts to locate the active analysis.
save(fullfile(project_dir,'current_config.mat'),'cfg','-v7.3');

parameter = ["project_dir";"analysis_label";"output_dir";"data_dir"; ...
    "n_subjects";"TR";"expected_timepoints"; ...
    "cv_scheme";"n_lambda";"voxel_block_size";"wm_smoothing"; ...
    "wm_smooth_fwhm_mm";"n_null";"subject_limit"];
value = [string(cfg.project_dir);string(cfg.analysis_label); ...
    string(cfg.output_dir);string(cfg.data_dir); ...
    string(numel(cfg.subject_list)); ...
    string(cfg.TR);string(cfg.expected_timepoints);string(cfg.cv_scheme); ...
    string(numel(cfg.lambda_grid));string(cfg.voxel_block_size); ...
    string(cfg.do_wm_mask_smoothing);string(cfg.wm_smooth_fwhm_mm); ...
    string(cfg.n_null);string(cfg.subject_limit)];
writetable(table(parameter,value), ...
    fullfile(cfg.result_dirs{1}, 'parameter_record.csv'));

fprintf('Configuration saved: %s\n', cfg.config_file);
fprintf('Subjects configured: %d\n', numel(cfg.subject_list));
fprintf('Main WM smoothing enabled: %d (FWHM %.1f mm)\n', ...
    cfg.do_wm_mask_smoothing, cfg.wm_smooth_fwhm_mm);

function name = local_stage_name(k)
names = { ...
 'config_project', ...
 'check_inputs_and_masks', ...
 'extract_timeseries', ...
 'build_model_matrices', ...
 'ridge_main_model', ...
 'write_r2_maps', ...
 'extract_beta_profiles', ...
 'write_roi_beta_maps', ...
 'network_contribution', ...
 'beta_similarity_gradient', ...
 'network_specificity', ...
 'fc_baseline', ...
 'network_level_ridge', ...
 'lagged_ridge', ...
 'null_model', ...
 'stability_analysis', ...
 'group_level_analysis', ...
 'behavior_prediction_optional', ...
 'summary_figures_tables'};
name = names{k+1};
end
