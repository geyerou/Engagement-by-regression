%% audit_existing_results
% Read-only audit of completed project results. New outputs are written to
% derivatives/result_audit_existing_results and original results are untouched.

clear; clc;
cfg = load_project_config();
out_dir = fullfile(cfg.output_dir,'result_audit_existing_results');
ensure_dir(out_dir);

wm = niftiread(cfg.wm_mask_final)>0;
wm_idx = find(wm);
n_sub = numel(cfg.subject_list);
n_vox = numel(wm_idx);

%% Main, network-level, and lagged model summaries
main_mean = nan(n_sub,1);
main_median = nan(n_sub,1);
main_positive = nan(n_sub,1);
main_corr = nan(n_sub,1);
net_mean = nan(n_sub,1);
lag_mean = nan(n_sub,1);
lag_delta_mean = nan(n_sub,1);
lag_delta_median = nan(n_sub,1);
lag_positive_delta = nan(n_sub,1);
best_lag_counts = zeros(n_sub,numel(cfg.lag_list_main));
lag_delta_group = zeros(n_sub,n_vox,'single');

for s = 1:n_sub
    sid = cfg.subject_list{s};
    M = load(fullfile(cfg.result_dirs{5},sprintf( ...
        'sub-%s_ridge_main_results.mat',sid)), ...
        'R2_cv_voxel','prediction_correlation');
    r2 = double(M.R2_cv_voxel(:));
    main_mean(s)=mean(r2); main_median(s)=median(r2);
    main_positive(s)=mean(r2>0);
    main_corr(s)=mean(double(M.prediction_correlation(:)));

    N = load(fullfile(cfg.result_dirs{13},sprintf( ...
        'sub-%s_network17_ridge_results.mat',sid)),'result');
    net_mean(s)=mean(double(N.result.R2(:)));

    L = load(fullfile(cfg.result_dirs{14},sprintf( ...
        'sub-%s_lagged_ridge_results.mat',sid)), ...
        'result','Delta_R2_lagged','best_lag','lags');
    lag_mean(s)=mean(double(L.result.R2(:)));
    d=double(L.Delta_R2_lagged(:));
    lag_delta_group(s,:)=single(d);
    lag_delta_mean(s)=mean(d);
    lag_delta_median(s)=median(d);
    lag_positive_delta(s)=mean(d>0);
    for j=1:numel(L.lags)
        best_lag_counts(s,j)=mean(double(L.best_lag(:))==L.lags(j));
    end
end

subject_id = string(cfg.subject_list(:));
model_subject = table(subject_id,main_mean,main_median,main_positive, ...
    main_corr,net_mean,lag_mean,lag_delta_mean,lag_delta_median, ...
    lag_positive_delta);
writetable(model_subject,fullfile(out_dir,'model_subject_audit.csv'));

lag_group_mean = mean(lag_delta_group,1)';
lag_group_consistency = mean(lag_delta_group>0,1)';
write_map_like(lag_group_mean,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_mean_lagged_minus_matched_zero_R2.nii.gz'));
write_map_like(lag_group_consistency,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_lagged_advantage_consistency.nii.gz'));

if license('test','statistics_toolbox')
    [~,p_lag,~,st_lag] = ttest(lag_delta_mean);
    lag_t = st_lag.tstat;
    p_main_vs_net = signrank(main_mean,net_mean);
    p_lag_vs_main = signrank(lag_mean,main_mean);
else
    lag_t = mean(lag_delta_mean)/(std(lag_delta_mean)/sqrt(n_sub));
    p_lag = NaN; p_main_vs_net=NaN; p_lag_vs_main=NaN;
end

%% Group R2 and significant extent
G = load(fullfile(cfg.result_dirs{17},'group_R2_statistics.mat'), ...
    'group_mean','group_median','q','consistency');
sig = G.q(:)<0.05 & G.group_mean(:)>0;

%% Network contribution and specificity
C = load(fullfile(cfg.result_dirs{17},'group_17net_contribution.mat'),'meanC');
one_net = dir(fullfile(cfg.result_dirs{9}, ...
    'sub-*_network_contribution_17net.mat'));
NN = load(fullfile(one_net(1).folder,one_net(1).name),'network_names');
network_names = string(NN.network_names(:));
mean_positive_delta = mean(max(double(C.meanC),0),1)';
mean_signed_delta = mean(double(C.meanC),1)';
positive_voxel_fraction = mean(double(C.meanC)>0,1)';
network_summary = table((1:17)',network_names,mean_signed_delta, ...
    mean_positive_delta,positive_voxel_fraction, ...
    'VariableNames',{'network_id','network_name','mean_signed_DeltaR2', ...
    'mean_positive_DeltaR2','positive_voxel_fraction'});
writetable(network_summary,fullfile(out_dir,'network_contribution_audit.csv'));

spec_map = niftiread(fullfile(cfg.result_dirs{17}, ...
    'group_mean_17net_specificity.nii.gz'));
spec = double(spec_map(wm_idx));

%% FC-beta correspondence and fingerprint stability
fc_beta_dir = fullfile(cfg.result_dirs{12},'s011b_fc_vs_beta');
fc_beta_signed_map = niftiread(fullfile(fc_beta_dir, ...
    'group_profile_FC_beta_signed_correlation.nii.gz'));
fc_beta_signed = double(fc_beta_signed_map(wm_idx));
stab_dir = fullfile(cfg.result_dirs{16},'s015a_beta_fc_stability');
beta_stab_map = niftiread(fullfile(stab_dir, ...
    'group_session_mean_beta_stability.nii.gz'));
fc_stab_map = niftiread(fullfile(stab_dir, ...
    'group_session_mean_FC_stability.nii.gz'));
beta_stab = double(beta_stab_map(wm_idx));
fc_stab = double(fc_stab_map(wm_idx));

%% Gradient results and cross-modal gradient agreement
grad_dir = fullfile(cfg.result_dirs{10},'s009a_brainspace_gradient');
GR = load(fullfile(grad_dir,'group_brainspace_gradient_results.mat'), ...
    'eigenvalues','gradients');
eig = double(GR.eigenvalues(:));
eig_fraction_first10 = eig(1:min(10,end))/sum(max(eig(1:min(10,end)),0));

gst_dir = fullfile(cfg.result_dirs{16},'s015b_gradient_stability');
session_beta_g1 = niftiread(fullfile(gst_dir, ...
    'session_beta_splitA_aligned_gradient-01.nii.gz'));
session_fc_g1 = niftiread(fullfile(gst_dir, ...
    'session_FC_splitA_aligned_gradient-01.nii.gz'));
beta_fc_gradient1_r = corr(double(session_beta_g1(wm_idx)), ...
    double(session_fc_g1(wm_idx)));

gradient_stability = readtable(fullfile(gst_dir, ...
    'beta_vs_FC_gradient_stability.csv'),'TextType','string');

%% Behavioral results
behavior = readtable(fullfile(cfg.result_dirs{18}, ...
    'behavior_prediction_summary.csv'),'TextType','string');
associations = readtable(fullfile(cfg.result_dirs{18}, ...
    'behavior_feature_partial_correlations.csv'),'TextType','string');

%% JHU anatomical localization
jhu_dir = fullfile(cfg.data_dir,'JHU');
jhu_file = fullfile(jhu_dir,'JHU-ICBM-labels-2mm.nii.gz');
jhu_xml = fullfile(jhu_dir,'JHU-ICBM-labels-2mm.xml');
jhu = niftiread(jhu_file);
if ~isequal(size(jhu),size(wm))
    error('JHU atlas and WM mask dimensions differ.');
end
jhu_info = niftiinfo(jhu_file);
voxel_mm = mean(double(jhu_info.PixelDimensions(1:3)));
xml_text = fileread(jhu_xml);
tok = regexp(xml_text, ...
    '<label index="(\d+)"[^>]*>(.*?)</label>','tokens');
n_label = numel(tok)-1;
label_id = zeros(n_label,1);
label_name = strings(n_label,1);
for k=2:numel(tok)
    label_id(k-1)=str2double(tok{k}{1});
    label_name(k-1)=strtrim(string(tok{k}{2}));
end

direct_label = double(jhu(wm_idx));
[distance_vox,nearest_linear] = bwdist(jhu>0);
nearest_distance_mm = double(distance_vox(wm_idx))*voxel_mm;
nearest_label = double(jhu(nearest_linear(wm_idx)));
nearest_label(nearest_distance_mm>6)=0;

jhu_rows = cell(n_label,18);
metric_names = ["group_mean_R2","lag_DeltaR2","specificity", ...
    "beta_stability","FC_stability","FC_beta_similarity"];
metric_data = [double(G.group_mean(:)),double(lag_group_mean),spec, ...
    beta_stab,fc_stab,fc_beta_signed];

for k=1:n_label
    id=label_id(k);
    direct = direct_label==id;
    nearest6 = nearest_label==id;
    dist_to_label = bwdist(jhu==id)*voxel_mm;
    dist_mask = double(dist_to_label(wm_idx));
    dilation4 = dist_mask<=4;
    min_distance = min(dist_mask);
    vals = nan(1,12);
    for m=1:size(metric_data,2)
        vals(m)=mean(metric_data(nearest6,m),'omitnan');
        vals(6+m)=mean(metric_data(dilation4,m),'omitnan');
    end
    jhu_rows(k,:) = [{id,label_name(k),sum(direct),sum(nearest6), ...
        sum(dilation4),min_distance},num2cell(vals(1:6)), ...
        num2cell(vals(7:12))];
end
jhu_names = {'label_id','label_name','direct_mask_voxels', ...
    'nearest_within_6mm_voxels','within_label_4mm_voxels', ...
    'minimum_mask_distance_mm'};
for m=1:numel(metric_names)
    jhu_names{end+1}=char("nearest6_"+metric_names(m)); %#ok<SAGROW>
end
for m=1:numel(metric_names)
    jhu_names{end+1}=char("dilated4_"+metric_names(m)); %#ok<SAGROW>
end
jhu_summary = cell2table(jhu_rows,'VariableNames',jhu_names);
writetable(jhu_summary,fullfile(out_dir,'JHU_tract_metric_audit.csv'));

coverage = table( ...
    sum(direct_label>0),mean(direct_label>0), ...
    sum(nearest_distance_mm<=2),mean(nearest_distance_mm<=2), ...
    sum(nearest_distance_mm<=4),mean(nearest_distance_mm<=4), ...
    sum(nearest_distance_mm<=6),mean(nearest_distance_mm<=6), ...
    median(nearest_distance_mm),prctile(nearest_distance_mm,95), ...
    'VariableNames',{'direct_labeled_voxels','direct_labeled_fraction', ...
    'within2mm_voxels','within2mm_fraction','within4mm_voxels', ...
    'within4mm_fraction','within6mm_voxels','within6mm_fraction', ...
    'median_distance_mm','p95_distance_mm'});
writetable(coverage,fullfile(out_dir,'JHU_mask_coverage.csv'));

nearest_vol = zeros(size(wm),'uint8');
nearest_vol(wm_idx)=uint8(nearest_label);
write_map_like(nearest_vol(wm_idx),wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'JHU_nearest_label_within_6mm.nii.gz'));
distance_vol = zeros(size(wm),'single');
distance_vol(wm_idx)=single(nearest_distance_mm);
write_map_like(distance_vol(wm_idx),wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'JHU_nearest_label_distance_mm.nii.gz'));

%% Compact global metric table
metric = [ ...
    "n_subjects";"n_wm_voxels";"mean_subject_main_R2"; ...
    "mean_subject_network17_R2";"mean_subject_lagged_R2"; ...
    "mean_subject_lagged_minus_matched_zero_R2"; ...
    "median_subject_lagged_minus_matched_zero_R2"; ...
    "fraction_subjects_mean_lag_advantage_positive"; ...
    "lag_advantage_t";"lag_advantage_p"; ...
    "main_vs_network_signrank_p";"lagged_vs_main_signrank_p"; ...
    "group_mean_R2";"group_positive_R2_fraction"; ...
    "group_FDR_significant_positive_fraction"; ...
    "mean_network_specificity";"mean_group_FC_beta_profile_r"; ...
    "mean_session_beta_fingerprint_stability"; ...
    "mean_session_FC_fingerprint_stability"; ...
    "beta_FC_gradient1_spatial_r_raw"; ...
    "beta_FC_gradient1_spatial_r_absolute"; ...
    "minimum_behavior_corrected_p"; ...
    "JHU_direct_coverage_fraction";"JHU_within6mm_coverage_fraction"];
value = [ ...
    n_sub;n_vox;mean(main_mean);mean(net_mean);mean(lag_mean); ...
    mean(lag_delta_mean);median(lag_delta_mean);mean(lag_delta_mean>0); ...
    lag_t;p_lag;p_main_vs_net;p_lag_vs_main; ...
    mean(double(G.group_mean));mean(double(G.group_mean)>0);mean(sig); ...
    mean(spec,'omitnan');mean(fc_beta_signed,'omitnan'); ...
    mean(beta_stab,'omitnan');mean(fc_stab,'omitnan'); ...
    beta_fc_gradient1_r;abs(beta_fc_gradient1_r); ...
    min([behavior.permutation_q_delta_r;behavior.permutation_q_delta_R2]); ...
    coverage.direct_labeled_fraction;coverage.within6mm_fraction];
global_metrics = table(metric,value);
writetable(global_metrics,fullfile(out_dir,'global_result_audit.csv'));

%% Figures
fig=figure('Color','w','Position',[40 40 1450 900],'Visible','off');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile;
model_labels = repmat(["400 ROI","17 network","Lagged"],n_sub,1);
boxchart(categorical(model_labels(:)),[main_mean;net_mean;lag_mean]);
ylabel('Mean voxelwise CV R^2'); title('Subject-level model performance');
grid on;
nexttile;
histogram(lag_delta_mean,20); xline(0,'--k');
xlabel('Subject mean lagged - matched zero-lag R^2');
title(sprintf('Lag gain: mean %.4g; %d/%d positive', ...
    mean(lag_delta_mean),sum(lag_delta_mean>0),n_sub));
nexttile;
bar(mean(best_lag_counts,1));
xticks(1:numel(cfg.lag_list_main)); xticklabels(cfg.lag_list_main);
xlabel('Best lag (TR)'); ylabel('Mean WM voxel fraction');
title('Lag-energy winning lag');
nexttile;
[~,ord]=sort(mean_positive_delta,'descend');
bar(mean_positive_delta(ord));
xticks(1:17); xticklabels(network_names(ord)); xtickangle(65);
ylabel('Mean positive Delta R^2'); title('Network ablation contribution');
nexttile;
scatter(fc_stab,beta_stab,3,fc_beta_signed,'filled');
hold on; lim=[min([xlim ylim]) max([xlim ylim])]; plot(lim,lim,'--k');
xlim(lim); ylim(lim); axis square; colorbar;
xlabel('FC fingerprint stability'); ylabel('Beta fingerprint stability');
title('Voxelwise reliability (color: FC-beta similarity)');
nexttile;
bar(100*eig_fraction_first10);
xlabel('Gradient'); ylabel('Relative eigenvalue (%)');
title(sprintf('Beta gradient spectrum; G1 beta-FC |r|=%.3f', ...
    abs(beta_fc_gradient1_r)));
exportgraphics(fig,fullfile(out_dir,'audit_overview.png'),'Resolution',200);
close(fig);

% Spatial overview uses representative axial slices through the common mask.
map_list = {double(G.group_mean(:)),'Group mean R^2'; ...
    double(lag_group_mean),'Lagged - zero-lag R^2'; ...
    spec,'Network specificity'; ...
    beta_stab,'Beta stability'; ...
    fc_stab,'FC stability'; ...
    fc_beta_signed,'FC-beta profile r'};
fig=figure('Color','w','Position',[40 40 1500 850],'Visible','off');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
z=46;
for m=1:size(map_list,1)
    vol=nan(size(wm),'single'); vol(wm_idx)=single(map_list{m,1});
    img=rot90(vol(:,:,z));
    nexttile;
    lo=prctile(map_list{m,1},2); hi=prctile(map_list{m,1},98);
    if lo==hi, hi=lo+eps; end
    imagesc(img,[lo hi]); axis image off; colorbar;
    colormap(gca,turbo(256)); title(map_list{m,2});
end
exportgraphics(fig,fullfile(out_dir,'audit_spatial_overview.png'), ...
    'Resolution',200);
close(fig);

% JHU tract display, sorted by direct/nearby group R2.
valid_jhu = jhu_summary.nearest_within_6mm_voxels>=20;
J = jhu_summary(valid_jhu,:);
[~,ord]=sort(J.nearest6_group_mean_R2,'descend');
J=J(ord,:);
nshow=min(15,height(J));
fig=figure('Color','w','Position',[40 40 1350 750],'Visible','off');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
barh(categorical(J.label_name(1:nshow)),J.nearest6_group_mean_R2(1:nshow));
set(gca,'YDir','reverse'); xlabel('Group mean R^2');
title('Top JHU tracts by prediction');
nexttile;
barh(categorical(J.label_name(1:nshow)), ...
    J.nearest6_lag_DeltaR2(1:nshow));
set(gca,'YDir','reverse'); xlabel('Lagged - zero-lag R^2');
title('Lag gain in the same tracts');
exportgraphics(fig,fullfile(out_dir,'JHU_tract_overview.png'), ...
    'Resolution',200);
close(fig);

save(fullfile(out_dir,'audit_workspace.mat'), ...
    'model_subject','global_metrics','network_summary','jhu_summary', ...
    'coverage','gradient_stability','behavior','associations', ...
    'lag_group_mean','lag_group_consistency','best_lag_counts', ...
    'eig','-v7.3');

fprintf('Result audit saved to: %s\n',out_dir);
