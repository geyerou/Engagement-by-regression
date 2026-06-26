%% s004b_refit_shared_lambda_betas
% Refit final standardized beta fingerprints using one common lambda.
%
% Prerequisite:
%   s004 must be complete for every configured subject.
%
% This stage does NOT replace the out-of-sample R2 from s004. It only creates
% a second set of full-data beta estimates with identical shrinkage strength
% across subjects, intended for cross-subject beta/profile comparisons.

clear; clc;
cfg = load_project_config();
input_dir = cfg.result_dirs{4};
ridge_dir = cfg.result_dirs{5};
out_dir = fullfile(ridge_dir, 'shared_lambda_refit');
ensure_dir(out_dir);

n_subjects = numel(cfg.subject_list);
if n_subjects < 2
    error(['Shared lambda requires multiple subjects. Complete the one-subject ' ...
        'debug run first, then set cfg.subject_limit = Inf, rerun s000, ' ...
        'complete s004 for all subjects, and run this script.']);
end

final_lambdas = zeros(n_subjects,1);
for s = 1:n_subjects
    subject_id = cfg.subject_list{s};
    result_file = fullfile(ridge_dir, ...
        sprintf('sub-%s_ridge_main_results.mat',subject_id));
    if ~exist(result_file,'file')
        error('Missing s004 result for subject %s: %s',subject_id,result_file);
    end
    D = load(result_file,'final_lambda');
    final_lambdas(s) = D.final_lambda;
end

switch lower(cfg.shared_lambda_method)
    case 'geometric_median'
        shared_lambda = exp(median(log(final_lambdas)));
    case 'geometric_mean'
        shared_lambda = exp(mean(log(final_lambdas)));
    otherwise
        error('Unknown shared_lambda_method: %s',cfg.shared_lambda_method);
end

lambda_table = table(string(cfg.subject_list(:)),final_lambdas, ...
    repmat(shared_lambda,n_subjects,1), ...
    'VariableNames',{'subject_id','subject_final_lambda','shared_lambda'});
writetable(lambda_table,fullfile(out_dir,'shared_lambda_record.csv'));
save(fullfile(out_dir,'shared_lambda_definition.mat'), ...
    'shared_lambda','final_lambdas','-v7.3');

fprintf('Shared lambda (%s) = %.8g\n', ...
    cfg.shared_lambda_method,shared_lambda);

% Subjects remain serial to control memory use.
for s = 1:n_subjects
    subject_id = cfg.subject_list{s};
    out_file = fullfile(out_dir, ...
        sprintf('sub-%s_shared_lambda_beta.mat',subject_id));
    if exist(out_file,'file')
        fprintf('Subject %s exists; skip.\n',subject_id);
        continue;
    end

    fprintf('Shared-lambda refit %s (%d/%d)\n', ...
        subject_id,s,n_subjects);
    model_input_file = fullfile(input_dir, ...
        sprintf('sub-%s_model_input.mat',subject_id));
    [X_runs,Y_runs,meta] = load_subject_runs(model_input_file);
    Xall = vertcat(X_runs{:});
    Yall = vertcat(Y_runs{:});

    [B_ridge_shared,muX_shared,sdX_shared,muY_shared,sdY_shared] = ...
        fit_final_beta_blocked(Xall,Yall,shared_lambda, ...
        cfg.voxel_block_size);

    wm_voxel_indices = meta.wm_voxel_indices;
    subject_final_lambda = final_lambdas(s);
    save(out_file,'subject_id','B_ridge_shared','shared_lambda', ...
        'subject_final_lambda','wm_voxel_indices','muX_shared', ...
        'sdX_shared','muY_shared','sdY_shared','-v7.3');

    clear X_runs Y_runs Xall Yall B_ridge_shared
end
