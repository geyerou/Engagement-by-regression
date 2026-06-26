%% s011b_compare_fc_vs_beta
% Compare traditional marginal FC fingerprints with conditional ridge beta
% fingerprints.
%
% Inputs
%   s006: shared-lambda Beta_profile, V x 400
%   s011: FC_profile, V x 400
%
% Outputs
%   - voxelwise signed-profile correlation
%   - voxelwise absolute-profile correlation
%   - 17-network composition correlation
%   - overlap of strongest absolute ROI weights
%   - subject and group NIfTI maps
%   - group mean FC/beta profiles and summary tables

clear; clc;
cfg = load_project_config();
labels = read_schaefer_labels(cfg);

beta_dir = cfg.result_dirs{7};
fc_dir = cfg.result_dirs{12};
out_dir = fullfile(fc_dir,'s011b_fc_vs_beta');
ensure_dir(out_dir);

n_subjects = numel(cfg.subject_list);
wm = niftiread(cfg.wm_mask_final)>0;
wm_idx = find(wm);
n_vox = numel(wm_idx);
n_roi = 400;
n_net = 17;
top_k = max(1,round(cfg.fc_beta_top_roi_fraction*n_roi));

signed_similarity = zeros(n_subjects,n_vox,'single');
absolute_similarity = zeros(n_subjects,n_vox,'single');
network_similarity = zeros(n_subjects,n_vox,'single');
top_roi_overlap = zeros(n_subjects,n_vox,'single');

group_fc_sum = zeros(n_vox,n_roi,'double');
group_beta_sum = zeros(n_vox,n_roi,'double');

summary_rows = cell(n_subjects,13);

for s = 1:n_subjects
    subject_id = cfg.subject_list{s};
    fprintf('FC-vs-beta subject %s (%d/%d)\n', ...
        subject_id,s,n_subjects);

    beta_file = fullfile(beta_dir, ...
        sprintf('sub-%s_beta_profile.mat',subject_id));
    fc_file = fullfile(fc_dir, ...
        sprintf('sub-%s_FC_profile.mat',subject_id));
    if ~exist(beta_file,'file')
        error('Missing s006 beta profile: %s',beta_file);
    end
    if ~exist(fc_file,'file')
        error('Missing s011 FC profile: %s',fc_file);
    end

    B = load(beta_file,'Beta_profile','wm_voxel_indices','beta_source');
    F = load(fc_file,'FC_profile','wm_voxel_indices');
    if ~isequal(B.wm_voxel_indices,F.wm_voxel_indices) || ...
            ~isequal(B.wm_voxel_indices,wm_idx)
        error('WM indices mismatch for subject %s.',subject_id);
    end
    if ~strcmp(B.beta_source,'shared_lambda')
        warning('Subject %s beta_source is %s, not shared_lambda.', ...
            subject_id,B.beta_source);
    end

    beta = single(B.Beta_profile);
    fc = single(F.FC_profile);
    if ~isequal(size(beta),[n_vox n_roi]) || ...
            ~isequal(size(fc),[n_vox n_roi])
        error('Unexpected FC/beta dimensions for subject %s.',subject_id);
    end
    if any(~isfinite(beta(:))) || any(~isfinite(fc(:)))
        error('FC or beta contains NaN/Inf for subject %s.',subject_id);
    end

    signed_r = row_correlation(fc,beta);
    absolute_r = row_correlation(abs(fc),abs(beta));

    fc_network = zeros(n_vox,n_net,'single');
    beta_network = zeros(n_vox,n_net,'single');
    for k = 1:n_net
        roi_mask = labels.network17_id==k;
        % Mean absolute weight corrects for unequal network parcel counts.
        fc_network(:,k) = mean(abs(fc(:,roi_mask)),2);
        beta_network(:,k) = mean(abs(beta(:,roi_mask)),2);
    end
    fc_network = fc_network ./ max(sum(fc_network,2),eps('single'));
    beta_network = beta_network ./ max(sum(beta_network,2),eps('single'));
    network_r = row_correlation(fc_network,beta_network);

    % Agreement between the strongest absolute FC and beta ROIs.
    [~,fc_top] = maxk(abs(fc),top_k,2);
    [~,beta_top] = maxk(abs(beta),top_k,2);
    overlap_count = zeros(n_vox,1,'single');
    for j = 1:top_k
        overlap_count = overlap_count + single(any( ...
            beta_top==fc_top(:,j),2));
    end
    overlap = overlap_count/top_k;

    signed_similarity(s,:) = signed_r;
    absolute_similarity(s,:) = absolute_r;
    network_similarity(s,:) = network_r;
    top_roi_overlap(s,:) = overlap;
    group_fc_sum = group_fc_sum + double(fc);
    group_beta_sum = group_beta_sum + double(beta);

    summary_rows(s,:) = {subject_id,B.beta_source, ...
        mean(signed_r,'omitnan'),median(signed_r,'omitnan'), ...
        mean(absolute_r,'omitnan'),median(absolute_r,'omitnan'), ...
        mean(network_r,'omitnan'),median(network_r,'omitnan'), ...
        mean(overlap,'omitnan'),median(overlap,'omitnan'), ...
        mean(signed_r>0.5),mean(network_r>0.5), ...
        mean(overlap>=0.5)};

    if cfg.write_subject_fc_beta_maps
        write_map_like(signed_r,wm_idx,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf( ...
            'sub-%s_FC_beta_signed_profile_correlation.nii.gz',subject_id)));
        write_map_like(absolute_r,wm_idx,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf( ...
            'sub-%s_FC_beta_absolute_profile_correlation.nii.gz',subject_id)));
        write_map_like(network_r,wm_idx,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf( ...
            'sub-%s_FC_beta_network17_correlation.nii.gz',subject_id)));
        write_map_like(overlap,wm_idx,cfg.wm_mask_final, ...
            fullfile(out_dir,sprintf( ...
            'sub-%s_FC_beta_topROI_overlap.nii.gz',subject_id)));
    end
    clear beta fc fc_network beta_network fc_top beta_top
end

subject_summary = cell2table(summary_rows,'VariableNames',{ ...
    'subject_id','beta_source', ...
    'mean_signed_profile_r','median_signed_profile_r', ...
    'mean_absolute_profile_r','median_absolute_profile_r', ...
    'mean_network17_r','median_network17_r', ...
    'mean_topROI_overlap','median_topROI_overlap', ...
    'fraction_signed_r_gt_05','fraction_network_r_gt_05', ...
    'fraction_topROI_overlap_ge_05'});
writetable(subject_summary,fullfile(out_dir, ...
    'FC_vs_beta_subject_summary.csv'));

group_mean_fc = single(group_fc_sum/n_subjects);
group_mean_beta = single(group_beta_sum/n_subjects);
group_signed_similarity = row_correlation(group_mean_fc,group_mean_beta);
group_absolute_similarity = row_correlation( ...
    abs(group_mean_fc),abs(group_mean_beta));

group_fc_network = zeros(n_vox,n_net,'single');
group_beta_network = zeros(n_vox,n_net,'single');
for k = 1:n_net
    roi_mask = labels.network17_id==k;
    group_fc_network(:,k) = mean(abs(group_mean_fc(:,roi_mask)),2);
    group_beta_network(:,k) = mean(abs(group_mean_beta(:,roi_mask)),2);
end
group_fc_network = group_fc_network ./ ...
    max(sum(group_fc_network,2),eps('single'));
group_beta_network = group_beta_network ./ ...
    max(sum(group_beta_network,2),eps('single'));
group_network_similarity = row_correlation( ...
    group_fc_network,group_beta_network);

[~,group_fc_top] = maxk(abs(group_mean_fc),top_k,2);
[~,group_beta_top] = maxk(abs(group_mean_beta),top_k,2);
group_overlap_count = zeros(n_vox,1,'single');
for j = 1:top_k
    group_overlap_count = group_overlap_count + single(any( ...
        group_beta_top==group_fc_top(:,j),2));
end
group_top_roi_overlap = group_overlap_count/top_k;

% Mean of subject-level correspondence maps answers a different question
% from correspondence between the two group-average profiles. Save both.
mean_subject_signed_similarity = mean(signed_similarity,1,'omitnan')';
mean_subject_absolute_similarity = mean(absolute_similarity,1,'omitnan')';
mean_subject_network_similarity = mean(network_similarity,1,'omitnan')';
mean_subject_top_roi_overlap = mean(top_roi_overlap,1,'omitnan')';

map_names = { ...
    'group_profile_FC_beta_signed_correlation',group_signed_similarity; ...
    'group_profile_FC_beta_absolute_correlation',group_absolute_similarity; ...
    'group_profile_FC_beta_network17_correlation',group_network_similarity; ...
    'group_profile_FC_beta_topROI_overlap',group_top_roi_overlap; ...
    'mean_subject_FC_beta_signed_correlation',mean_subject_signed_similarity; ...
    'mean_subject_FC_beta_absolute_correlation',mean_subject_absolute_similarity; ...
    'mean_subject_FC_beta_network17_correlation',mean_subject_network_similarity; ...
    'mean_subject_FC_beta_topROI_overlap',mean_subject_top_roi_overlap};
for i = 1:size(map_names,1)
    write_map_like(map_names{i,2},wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,[map_names{i,1} '.nii.gz']));
end

save(fullfile(out_dir,'FC_vs_beta_group_results.mat'), ...
    'group_mean_fc','group_mean_beta','group_fc_network', ...
    'group_beta_network','group_signed_similarity', ...
    'group_absolute_similarity','group_network_similarity', ...
    'group_top_roi_overlap','signed_similarity','absolute_similarity', ...
    'network_similarity','top_roi_overlap','wm_idx','top_k', ...
    'labels','subject_summary','-v7.3');

group_metric = ["signed profile correlation"; ...
    "absolute profile correlation"; ...
    "17-network composition correlation"; ...
    "top ROI overlap"];
group_profile_value = [mean(group_signed_similarity,'omitnan'); ...
    mean(group_absolute_similarity,'omitnan'); ...
    mean(group_network_similarity,'omitnan'); ...
    mean(group_top_roi_overlap,'omitnan')];
mean_subject_value = [mean(mean_subject_signed_similarity,'omitnan'); ...
    mean(mean_subject_absolute_similarity,'omitnan'); ...
    mean(mean_subject_network_similarity,'omitnan'); ...
    mean(mean_subject_top_roi_overlap,'omitnan')];
group_summary = table(group_metric,group_profile_value,mean_subject_value);
writetable(group_summary,fullfile(out_dir,'FC_vs_beta_group_summary.csv'));

fig = figure('Color','w','Position',[50 50 1300 750],'Visible','off');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; histogram(group_signed_similarity,100);
xlabel('Signed FC-beta profile correlation'); ylabel('WM voxels');
title(sprintf('Group profiles, mean r = %.3f', ...
    mean(group_signed_similarity,'omitnan')));
nexttile; histogram(group_absolute_similarity,100);
xlabel('Absolute-pattern correlation'); ylabel('WM voxels');
title(sprintf('Group profiles, mean r = %.3f', ...
    mean(group_absolute_similarity,'omitnan')));
nexttile; histogram(group_network_similarity,100);
xlabel('17-network composition correlation'); ylabel('WM voxels');
title(sprintf('Group profiles, mean r = %.3f', ...
    mean(group_network_similarity,'omitnan')));
nexttile; histogram(group_top_roi_overlap,0:0.05:1);
xlabel(sprintf('Top-%d ROI overlap',top_k)); ylabel('WM voxels');
title(sprintf('Group profiles, mean overlap = %.3f', ...
    mean(group_top_roi_overlap,'omitnan')));
exportgraphics(fig,fullfile(out_dir,'FC_vs_beta_group_distributions.png'), ...
    'Resolution',200);
close(fig);

fprintf('FC-vs-beta comparison saved to: %s\n',out_dir);
