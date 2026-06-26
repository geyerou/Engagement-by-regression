%% s015a_compare_beta_fc_stability
% Compare split-half reliability of ridge beta and marginal FC fingerprints.
%
% Two complementary splits are used:
%   session  : runs [1 2] versus [3 4]
%   balanced : runs [1 4] versus [2 3]
%
% Outputs also include group-average split profiles used by s015b.

clear; clc;
cfg = load_project_config();
input_dir = cfg.result_dirs{4};
ridge_dir = cfg.result_dirs{5};
out_dir = fullfile(cfg.result_dirs{16},'s015a_beta_fc_stability');
ensure_dir(out_dir);

n_subjects = numel(cfg.subject_list);
wm = niftiread(cfg.wm_mask_final)>0;
wm_idx = find(wm);
n_vox = numel(wm_idx);
n_roi = 400;

split_names = string(cfg.stability_split_names);
split_runs = {
    {[1 2],[3 4]}, ...
    {[1 4],[2 3]}};
if numel(split_names) ~= numel(split_runs)
    error('stability_split_names does not match the implemented splits.');
end

if cfg.stability_use_shared_lambda
    shared_file = fullfile(ridge_dir,'shared_lambda_refit', ...
        'shared_lambda_definition.mat');
    if ~exist(shared_file,'file')
        error('Shared lambda definition missing. Run s004b first.');
    end
    L = load(shared_file,'shared_lambda');
    stability_lambda = L.shared_lambda;
    lambda_source = 'shared_lambda';
else
    stability_lambda = NaN;
    lambda_source = 'subject_final_lambda';
end

summary_rows = cell(n_subjects*numel(split_names),12);
row_counter = 0;

for q = 1:numel(split_names)
    split_name = split_names(q);
    runs_a = split_runs{q}{1};
    runs_b = split_runs{q}{2};

    beta_stability = zeros(n_subjects,n_vox,'single');
    fc_stability = zeros(n_subjects,n_vox,'single');
    beta_minus_fc = zeros(n_subjects,n_vox,'single');

    group_beta_a_sum = zeros(n_vox,n_roi,'double');
    group_beta_b_sum = zeros(n_vox,n_roi,'double');
    group_fc_a_sum = zeros(n_vox,n_roi,'double');
    group_fc_b_sum = zeros(n_vox,n_roi,'double');

    for s = 1:n_subjects
        subject_id = cfg.subject_list{s};
        fprintf('%s split subject %s (%d/%d)\n', ...
            split_name,subject_id,s,n_subjects);
        [X_runs,Y_runs,meta] = load_subject_runs(fullfile(input_dir, ...
            sprintf('sub-%s_model_input.mat',subject_id)));
        if ~isequal(meta.wm_voxel_indices,wm_idx)
            error('WM indices mismatch for subject %s.',subject_id);
        end

        if cfg.stability_use_shared_lambda
            lambda = stability_lambda;
        else
            M = load(fullfile(ridge_dir,sprintf( ...
                'sub-%s_ridge_main_results.mat',subject_id)), ...
                'final_lambda');
            lambda = M.final_lambda;
        end

        Xa = vertcat(X_runs{runs_a});
        Ya = vertcat(Y_runs{runs_a});
        Xb = vertcat(X_runs{runs_b});
        Yb = vertcat(Y_runs{runs_b});

        [Ba,~,~,~,~] = fit_final_beta_blocked( ...
            Xa,Ya,lambda,cfg.voxel_block_size);
        [Bb,~,~,~,~] = fit_final_beta_blocked( ...
            Xb,Yb,lambda,cfg.voxel_block_size);

        FCa = compute_fc_profile(Xa,Ya,cfg.voxel_block_size);
        FCb = compute_fc_profile(Xb,Yb,cfg.voxel_block_size);

        beta_r = row_correlation(Ba',Bb');
        fc_r = row_correlation(FCa,FCb);
        difference = beta_r-fc_r;

        beta_stability(s,:) = beta_r;
        fc_stability(s,:) = fc_r;
        beta_minus_fc(s,:) = difference;

        group_beta_a_sum = group_beta_a_sum+double(Ba');
        group_beta_b_sum = group_beta_b_sum+double(Bb');
        group_fc_a_sum = group_fc_a_sum+double(FCa);
        group_fc_b_sum = group_fc_b_sum+double(FCb);

        row_counter = row_counter+1;
        summary_rows(row_counter,:) = {subject_id,split_name,lambda, ...
            mean(beta_r,'omitnan'),median(beta_r,'omitnan'), ...
            mean(fc_r,'omitnan'),median(fc_r,'omitnan'), ...
            mean(difference,'omitnan'),median(difference,'omitnan'), ...
            mean(beta_r>fc_r),mean(beta_r>0.5),mean(fc_r>0.5)};

        if cfg.write_subject_stability_maps
            write_map_like(beta_r,wm_idx,cfg.wm_mask_final, ...
                fullfile(out_dir,sprintf( ...
                'sub-%s_%s_beta_stability.nii.gz',subject_id,split_name)));
            write_map_like(fc_r,wm_idx,cfg.wm_mask_final, ...
                fullfile(out_dir,sprintf( ...
                'sub-%s_%s_FC_stability.nii.gz',subject_id,split_name)));
            write_map_like(difference,wm_idx,cfg.wm_mask_final, ...
                fullfile(out_dir,sprintf( ...
                'sub-%s_%s_beta_minus_FC_stability.nii.gz', ...
                subject_id,split_name)));
        end
        clear X_runs Y_runs Xa Xb Ya Yb Ba Bb FCa FCb
    end

    group_beta_split_a = single(group_beta_a_sum/n_subjects);
    group_beta_split_b = single(group_beta_b_sum/n_subjects);
    group_fc_split_a = single(group_fc_a_sum/n_subjects);
    group_fc_split_b = single(group_fc_b_sum/n_subjects);

    mean_beta_stability = mean(beta_stability,1,'omitnan')';
    mean_fc_stability = mean(fc_stability,1,'omitnan')';
    mean_beta_minus_fc = mean(beta_minus_fc,1,'omitnan')';
    write_map_like(mean_beta_stability,wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf( ...
        'group_%s_mean_beta_stability.nii.gz',split_name)));
    write_map_like(mean_fc_stability,wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf( ...
        'group_%s_mean_FC_stability.nii.gz',split_name)));
    write_map_like(mean_beta_minus_fc,wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf( ...
        'group_%s_mean_beta_minus_FC_stability.nii.gz',split_name)));

    save(fullfile(out_dir,sprintf( ...
        'group_%s_split_profiles.mat',split_name)), ...
        'group_beta_split_a','group_beta_split_b', ...
        'group_fc_split_a','group_fc_split_b','runs_a','runs_b', ...
        'stability_lambda','lambda_source','wm_idx','-v7.3');
    save(fullfile(out_dir,sprintf( ...
        'group_%s_subject_stability.mat',split_name)), ...
        'beta_stability','fc_stability','beta_minus_fc','wm_idx','-v7.3');
end

subject_summary = cell2table(summary_rows,'VariableNames',{ ...
    'subject_id','split_name','lambda', ...
    'mean_beta_stability','median_beta_stability', ...
    'mean_FC_stability','median_FC_stability', ...
    'mean_beta_minus_FC','median_beta_minus_FC', ...
    'fraction_voxels_beta_gt_FC','fraction_beta_r_gt_05', ...
    'fraction_FC_r_gt_05'});
writetable(subject_summary,fullfile(out_dir, ...
    'beta_vs_FC_stability_subject_summary.csv'));

group_rows = cell(numel(split_names),8);
for q = 1:numel(split_names)
    split_name = split_names(q);
    D = load(fullfile(out_dir,sprintf( ...
        'group_%s_subject_stability.mat',split_name)));
    group_rows(q,:) = {split_name, ...
        mean(D.beta_stability,'all','omitnan'), ...
        median(D.beta_stability,'all','omitnan'), ...
        mean(D.fc_stability,'all','omitnan'), ...
        median(D.fc_stability,'all','omitnan'), ...
        mean(D.beta_minus_fc,'all','omitnan'), ...
        mean(D.beta_stability>D.fc_stability,'all'), ...
        mean(D.beta_minus_fc>0.1,'all')};
end
group_summary = cell2table(group_rows,'VariableNames',{ ...
    'split_name','mean_beta_stability','median_beta_stability', ...
    'mean_FC_stability','median_FC_stability','mean_beta_minus_FC', ...
    'fraction_beta_gt_FC','fraction_beta_advantage_gt_01'});
writetable(group_summary,fullfile(out_dir, ...
    'beta_vs_FC_stability_group_summary.csv'));

fig = figure('Color','w','Position',[50 50 1200 500],'Visible','off');
tiledlayout(1,numel(split_names),'TileSpacing','compact');
for q = 1:numel(split_names)
    nexttile;
    rows = subject_summary.split_name==split_names(q);
    scatter(subject_summary.mean_FC_stability(rows), ...
        subject_summary.mean_beta_stability(rows),30,'filled');
    hold on;
    lim = [min([xlim ylim]) max([xlim ylim])];
    plot(lim,lim,'--k'); xlim(lim); ylim(lim); axis square;
    xlabel('Mean FC stability'); ylabel('Mean beta stability');
    title(sprintf('%s split',split_names(q)));
    grid on;
end
exportgraphics(fig,fullfile(out_dir, ...
    'beta_vs_FC_subject_stability.png'),'Resolution',200);
close(fig);

fprintf('Beta-vs-FC stability results saved to: %s\n',out_dir);
