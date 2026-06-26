%% s009b_beta_kmeans_parcellation
% Discrete WM parcellation directly from group mean shared-lambda beta
% fingerprints. This stage is independent of s009a gradients.
%
% PCA is used because the 400 ridge coefficients are highly correlated.
% Components are retained by explained variance rather than a fixed
% arbitrary dimension.

clear; clc;
cfg = load_project_config();
if ~cfg.do_beta_parcellation
    fprintf('Beta parcellation disabled in config.\n');
    return;
end
if ~license('test','statistics_toolbox')
    error('This stage requires Statistics and Machine Learning Toolbox.');
end

stage_dir = cfg.result_dirs{10};
out_dir = fullfile(stage_dir,'s009b_beta_kmeans_parcellation');
ensure_dir(out_dir);
rng(cfg.parcellation_random_state,'twister');

[profile,wm_idx] = group_mean_beta_profile(cfg);
profile = double(profile);
if any(~isfinite(profile(:)))
    error('Group beta profile contains NaN or Inf.');
end

% Cluster fingerprint shape rather than total beta magnitude.
profile = profile - mean(profile,2);
row_norm = sqrt(sum(profile.^2,2));
if any(row_norm <= eps)
    error('%d WM voxels have a zero centered beta fingerprint.', ...
        sum(row_norm<=eps));
end
profile = profile ./ row_norm;

max_pc = min([cfg.parcellation_pca_max_components, ...
    size(profile,2),size(profile,1)-1]);
[coeff,score,latent,~,explained,mu] = pca(profile, ...
    'NumComponents',max_pc,'Centered',true);
% Some MATLAB releases return a full-length explained vector even when
% score/coeff are truncated by NumComponents. Only components that actually
% exist in score can be selected for clustering.
n_available_pc = size(score,2);
explained_available = explained(1:n_available_pc);
cum_explained = cumsum(explained_available);
n_pc = find(cum_explained >= ...
    cfg.parcellation_pca_variance_percent,1,'first');
if isempty(n_pc)
    n_pc = n_available_pc;
    warning(['The first %d available PCs explain %.2f%% variance, below the ' ...
        'requested %.2f%%. All available PCs will be used.'], ...
        n_available_pc,cum_explained(end), ...
        cfg.parcellation_pca_variance_percent);
end
X = score(:,1:n_pc);
fprintf('PCA retained %d components (%.2f%% variance).\n', ...
    n_pc,cum_explained(n_pc));

parcellations = struct();
summary_rows = cell(numel(cfg.parcellation_k),6);
mask = niftiread(cfg.wm_mask_final)>0;
for q = 1:numel(cfg.parcellation_k)
    k = cfg.parcellation_k(q);
    fprintf('Running k-means with k=%d, replicates=%d.\n', ...
        k,cfg.parcellation_replicates);
    [cluster_id,centroids,sumd] = kmeans(X,k, ...
        'Replicates',cfg.parcellation_replicates, ...
        'MaxIter',cfg.parcellation_max_iter, ...
        'Start','plus','EmptyAction','singleton', ...
        'Options',statset('UseParallel',false));
    counts = accumarray(cluster_id,1,[k 1]);

    sample_n = min(cfg.parcellation_silhouette_sample,size(X,1));
    sample_idx = round(linspace(1,size(X,1),sample_n));
    sil = silhouette(X(sample_idx,:),cluster_id(sample_idx),'sqeuclidean');
    mean_silhouette = mean(sil);
    field_name = sprintf('k%d',k);
    parcellations.(field_name) = struct( ...
        'cluster_id',int16(cluster_id), ...
        'centroids',single(centroids), ...
        'within_cluster_sum',sumd, ...
        'cluster_sizes',counts, ...
        'sample_indices',sample_idx, ...
        'sample_silhouette',single(sil), ...
        'mean_silhouette',mean_silhouette);

    write_map_like(cluster_id,wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('group_beta_kmeans_k%02d.nii.gz',k)));
    summary_rows(q,:) = {k,n_pc,cum_explained(n_pc), ...
        mean_silhouette,min(counts),max(counts)};

    fig = figure('Color','w','Position',[50 50 1300 700],'Visible','off');
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    nexttile; bar(counts);
    xlabel('Cluster'); ylabel('WM voxels');
    title(sprintf('Cluster sizes, k=%d',k));
    nexttile; histogram(sil,50);
    xlabel('Silhouette'); ylabel('Sample voxels');
    title(sprintf('Mean silhouette = %.3f',mean_silhouette));
    z_slices = [38 46 54];
    vol = zeros(size(mask),'single'); vol(wm_idx)=cluster_id;
    for j=1:3
        nexttile;
        img=rot90(vol(:,:,z_slices(j)));
        m=rot90(mask(:,:,z_slices(j))); img(~m)=NaN;
        imagesc(img,[1 k]); axis image off;
        colormap(gca,lines(k)); colorbar;
        title(sprintf('k=%d, axial z=%d',k,z_slices(j)));
    end
    exportgraphics(fig,fullfile(out_dir,sprintf( ...
        'beta_kmeans_k%02d_summary.png',k)),'Resolution',200);
    close(fig);
end

summary = cell2table(summary_rows,'VariableNames', ...
    {'k','n_pca_components','variance_explained_percent', ...
    'mean_sample_silhouette','smallest_cluster','largest_cluster'});
writetable(summary,fullfile(out_dir,'beta_kmeans_summary.csv'));
save(fullfile(out_dir,'group_beta_kmeans_results.mat'), ...
    'parcellations','coeff','latent','explained','cum_explained','mu', ...
    'n_pc','wm_idx','summary','-v7.3');

fig = figure('Color','w','Visible','off');
plot(1:numel(cum_explained),cum_explained,'LineWidth',1.5);
hold on; yline(cfg.parcellation_pca_variance_percent,'--');
xline(n_pc,'--'); grid on;
xlabel('PCA component'); ylabel('Cumulative explained variance (%)');
title('Beta fingerprint PCA');
exportgraphics(fig,fullfile(out_dir,'beta_pca_variance.png'), ...
    'Resolution',200);
close(fig);
fprintf('Beta k-means parcellations saved to: %s\n',out_dir);
