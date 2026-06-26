%% s017_behavior_prediction
% Complete cognition analysis using ridge-derived WM features.
%
% Primary outcomes:
%   CogTotalComp_AgeAdj, CogFluidComp_AgeAdj,
%   CogCrystalComp_AgeAdj, PMAT24_A_CR
%
% Confounds:
%   age-group dummy variables, sex, and mean motion across four runs.
%
% Prediction:
%   repeated nested cross-validation compares a confound-only model against
%   confounds plus WM imaging features. All residualization, standardization,
%   and lambda selection occur within training folds.
%
% Inference:
%   Freedman-Lane residual permutation tests the incremental predictive value
%   of imaging features beyond confounds.

clear; clc;
cfg = load_project_config();
rng(cfg.behavior_random_state,'twister');

ridge_dir = cfg.result_dirs{5};
net_dir = cfg.result_dirs{9};
spec_dir = cfg.result_dirs{11};
out_dir = cfg.result_dirs{18};
ensure_dir(out_dir);

n = numel(cfg.subject_list);
subject_id = string(cfg.subject_list(:));

% -------------------------------------------------------------------------
% 1. Build imaging feature matrix.
% -------------------------------------------------------------------------
mean_R2 = zeros(n,1);
median_R2 = zeros(n,1);
positive_R2_fraction = zeros(n,1);
R2_p90 = zeros(n,1);
mean_specificity_17 = zeros(n,1);
network_contribution_17 = zeros(n,17);

for s = 1:n
    sid = subject_id(s);
    ridge_file = fullfile(ridge_dir,sprintf( ...
        'sub-%s_ridge_main_results.mat',sid));
    network_file = fullfile(net_dir,sprintf( ...
        'sub-%s_network_contribution_17net.mat',sid));
    specificity_file = fullfile(spec_dir,sprintf( ...
        'sub-%s_network_specificity_17net.mat',sid));
    if ~exist(ridge_file,'file') || ~exist(network_file,'file') || ...
            ~exist(specificity_file,'file')
        error(['s017 requires completed s004, s008, and s010 outputs. ' ...
            'Missing files for subject %s.'],sid);
    end

    R = load(ridge_file,'R2_cv_voxel');
    r2 = double(R.R2_cv_voxel(:));
    mean_R2(s) = mean(r2);
    median_R2(s) = median(r2);
    positive_R2_fraction(s) = mean(r2>0);
    R2_p90(s) = prctile(r2,90);

    C = load(network_file,'Delta_R2');
    network_contribution_17(s,:) = mean(max(double(C.Delta_R2),0),1);

    S = load(specificity_file,'Specificity');
    mean_specificity_17(s) = mean(double(S.Specificity),'omitnan');
end

feature_table = table(subject_id,mean_R2,median_R2, ...
    positive_R2_fraction,R2_p90,mean_specificity_17);
for k = 1:17
    feature_table.(sprintf('network17_contribution_%02d',k)) = ...
        network_contribution_17(:,k);
end
feature_names = string(feature_table.Properties.VariableNames(2:end));
X = table2array(feature_table(:,2:end));

% -------------------------------------------------------------------------
% 2. Join behavioral and confound data.
% -------------------------------------------------------------------------
B = readtable(cfg.behavior_table,'TextType','string');
B.Subject = string(B.Subject);
M = readtable(cfg.motion_table,'TextType','string');
M.Subject = string(M.Subject);

[found_b,loc_b] = ismember(subject_id,B.Subject);
[found_m,loc_m] = ismember(subject_id,M.Subject);
if ~all(found_b) || ~all(found_m)
    error('Behavior or motion table does not cover all configured subjects.');
end
B = B(loc_b,:);
M = M(loc_m,:);

age_group = categorical(B.Age);
sex = categorical(B.Gender);
age_levels = categories(age_group);
sex_levels = categories(sex);

% Reference-coded confounds: intercept + age dummies + sex dummies + motion.
C = ones(n,1);
confound_names = "intercept";
for k = 2:numel(age_levels)
    C(:,end+1) = double(age_group==age_levels{k}); %#ok<SAGROW>
    confound_names(end+1) = "age_"+string(age_levels{k}); %#ok<SAGROW>
end
for k = 2:numel(sex_levels)
    C(:,end+1) = double(sex==sex_levels{k}); %#ok<SAGROW>
    confound_names(end+1) = "sex_"+string(sex_levels{k}); %#ok<SAGROW>
end
motion_columns = {'rfMRI_REST1_LR','rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR','rfMRI_REST2_RL'};
motion_matrix = zeros(n,4);
for k = 1:4
    motion_matrix(:,k) = double(M.(motion_columns{k}));
end
mean_motion = mean(motion_matrix,2);
C(:,end+1) = mean_motion;
confound_names(end+1) = "mean_motion";

% Save the complete analysis table.
analysis_table = feature_table;
analysis_table.Age = string(B.Age);
analysis_table.Gender = string(B.Gender);
analysis_table.mean_motion = mean_motion;
for k = 1:numel(cfg.behavior_outcomes)
    outcome_name = cfg.behavior_outcomes{k};
    analysis_table.(outcome_name) = double(B.(outcome_name));
end
writetable(analysis_table,fullfile(out_dir, ...
    'behavior_prediction_analysis_table.csv'));

% Fixed repeated stratified outer folds are reused for observed and
% permutation analyses.
strata = string(B.Age)+"_"+string(B.Gender);
outer_folds = make_repeated_stratified_folds(strata, ...
    cfg.behavior_outer_folds,cfg.behavior_outer_repeats, ...
    cfg.behavior_random_state);

% -------------------------------------------------------------------------
% 3. Descriptive confound-adjusted feature/outcome associations.
% -------------------------------------------------------------------------
association_rows = cell(numel(cfg.behavior_outcomes)*size(X,2),7);
association_counter = 0;

% -------------------------------------------------------------------------
% 4. Nested prediction and Freedman-Lane permutation inference.
% -------------------------------------------------------------------------
prediction_rows = cell(numel(cfg.behavior_outcomes),15);
prediction_store = struct();

for o = 1:numel(cfg.behavior_outcomes)
    outcome_name = cfg.behavior_outcomes{o};
    y = double(B.(outcome_name));
    valid = isfinite(y) & all(isfinite(X),2) & all(isfinite(C),2);
    if nnz(valid) < 50
        error('Outcome %s has only %d complete observations.', ...
            outcome_name,nnz(valid));
    end

    yv = y(valid);
    Xv = X(valid,:);
    Cv = C(valid,:);
    fold_v = outer_folds(valid,:);
    sid_v = subject_id(valid);

    % Descriptive partial correlations. These are not predictive metrics.
    cy = Cv\yv;
    y_res = yv-Cv*cy;
    for f = 1:size(Xv,2)
        bx = Cv\Xv(:,f);
        x_res = Xv(:,f)-Cv*bx;
        pearson_r = corr(x_res,y_res,'Type','Pearson');
        spearman_r = corr(x_res,y_res,'Type','Spearman');
        [~,p_pearson] = corr(x_res,y_res,'Type','Pearson');
        [~,p_spearman] = corr(x_res,y_res,'Type','Spearman');
        association_counter = association_counter+1;
        association_rows(association_counter,:) = {outcome_name, ...
            feature_names(f),pearson_r,p_pearson,spearman_r,p_spearman, ...
            nnz(valid)};
    end

    fprintf('Nested behavior prediction: %s (N=%d)\n', ...
        outcome_name,nnz(valid));
    observed = nested_behavior_prediction(Xv,Cv,yv,fold_v, ...
        cfg.behavior_inner_folds,cfg.behavior_lambda_grid, ...
        cfg.behavior_random_state+o);

    % Freedman-Lane: retain confound fit, permute confound residuals.
    confound_fit_all = Cv*(Cv\yv);
    confound_residual_all = yv-confound_fit_all;
    null_delta_R2 = zeros(cfg.behavior_n_permutations,1);
    null_delta_r = zeros(cfg.behavior_n_permutations,1);
    for p = 1:cfg.behavior_n_permutations
        y_perm = confound_fit_all + ...
            confound_residual_all(randperm(numel(yv)));
        % Use the fold-specific penalties selected by the observed nested
        % analysis. This keeps every permuted prediction out-of-sample while
        % making 1000 permutations computationally feasible.
        perm_result = fixed_lambda_behavior_prediction( ...
            Xv,Cv,y_perm,fold_v,observed.selected_lambda);
        null_delta_R2(p) = perm_result.delta_R2;
        null_delta_r(p) = perm_result.delta_r;
        if mod(p,100)==0
            fprintf('  permutation %d/%d\n',p,cfg.behavior_n_permutations);
        end
    end
    p_delta_R2 = (1+sum(null_delta_R2>=observed.delta_R2))/ ...
        (1+cfg.behavior_n_permutations);
    p_delta_r = (1+sum(null_delta_r>=observed.delta_r))/ ...
        (1+cfg.behavior_n_permutations);

    prediction_rows(o,:) = {outcome_name,nnz(valid), ...
        observed.full_r,observed.full_R2,observed.full_MAE, ...
        observed.confound_r,observed.confound_R2,observed.confound_MAE, ...
        observed.delta_r,observed.delta_R2,p_delta_r,p_delta_R2, ...
        median(observed.selected_lambda), ...
        prctile(observed.selected_lambda,25), ...
        prctile(observed.selected_lambda,75)};

    field_name = matlab.lang.makeValidName(outcome_name);
    prediction_store.(field_name) = struct( ...
        'subject_id',sid_v,'observed',observed, ...
        'null_delta_R2',null_delta_R2,'null_delta_r',null_delta_r, ...
        'p_delta_R2',p_delta_R2,'p_delta_r',p_delta_r);
end

association_table = cell2table(association_rows,'VariableNames',{ ...
    'outcome','feature','partial_pearson_r','partial_pearson_p', ...
    'partial_spearman_r','partial_spearman_p','N'});
association_table.partial_pearson_q = fdr_bh( ...
    association_table.partial_pearson_p);
association_table.partial_spearman_q = fdr_bh( ...
    association_table.partial_spearman_p);
writetable(association_table,fullfile(out_dir, ...
    'behavior_feature_partial_correlations.csv'));

prediction_summary = cell2table(prediction_rows,'VariableNames',{ ...
    'outcome','N','full_model_r','full_model_R2','full_model_MAE', ...
    'confound_model_r','confound_model_R2','confound_model_MAE', ...
    'delta_r','delta_R2','permutation_p_delta_r', ...
    'permutation_p_delta_R2','median_selected_lambda', ...
    'lambda_p25','lambda_p75'});
prediction_summary.permutation_q_delta_r = fdr_bh( ...
    prediction_summary.permutation_p_delta_r);
prediction_summary.permutation_q_delta_R2 = fdr_bh( ...
    prediction_summary.permutation_p_delta_R2);
writetable(prediction_summary,fullfile(out_dir, ...
    'behavior_prediction_summary.csv'));

save(fullfile(out_dir,'behavior_prediction_results.mat'), ...
    'prediction_store','prediction_summary','association_table', ...
    'analysis_table','feature_names','confound_names','outer_folds','-v7.3');

fig = figure('Color','w','Position',[50 50 1200 550],'Visible','off');
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
bar(categorical(prediction_summary.outcome), ...
    [prediction_summary.confound_model_R2 prediction_summary.full_model_R2]);
yline(0,'--'); ylabel('Cross-validated R^2');
legend({'Confounds only','Confounds + imaging'},'Location','best');
title('Nested cross-validated cognition prediction'); grid on;
nexttile;
bar(categorical(prediction_summary.outcome),prediction_summary.delta_R2);
yline(0,'--'); ylabel('\DeltaR^2 from imaging');
title('Incremental imaging prediction'); grid on;
exportgraphics(fig,fullfile(out_dir, ...
    'behavior_prediction_summary.png'),'Resolution',200);
close(fig);

fprintf('Behavior analysis saved to: %s\n',out_dir);
