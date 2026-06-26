%% s009a_brainspace_gradient
% Standard BrainSpace gradient analysis of group mean shared-lambda beta
% fingerprints. Each WM voxel is a seed and its 400 ridge coefficients are
% the feature profile.
%
% BrainSpace constructs a dense 32445 x 32445 affinity matrix. Run this
% stage alone, without another memory-intensive MATLAB job.

clear; clc;
cfg = load_project_config();
if ~cfg.do_gradient
    fprintf('Gradient analysis disabled in config.\n');
    return;
end

if ~exist(cfg.brainspace_path,'dir')
    error('BrainSpace path not found: %s',cfg.brainspace_path);
end
addpath(genpath(cfg.brainspace_path));
if isempty(which('GradientMaps'))
    error('BrainSpace GradientMaps is unavailable after adding the path.');
end

stage_dir = cfg.result_dirs{10};
out_dir = fullfile(stage_dir,'s009a_brainspace_gradient');
ensure_dir(out_dir);

[profile,wm_idx] = group_mean_beta_profile(cfg);
profile = double(profile);
if any(~isfinite(profile(:)))
    error('Group beta profile contains NaN or Inf.');
end
row_norm = sqrt(sum(profile.^2,2));
if any(row_norm <= eps)
    error('%d WM voxels have a zero beta fingerprint.',sum(row_norm<=eps));
end

n_vox = size(profile,1);
affinity_gib = n_vox^2 * 8 / 1024^3;
fprintf('BrainSpace input: %d WM voxels x %d beta features.\n', ...
    size(profile,1),size(profile,2));
fprintf('One dense double affinity matrix requires approximately %.2f GiB.\n', ...
    affinity_gib);
fprintf('Peak memory will be higher because intermediate matrices are created.\n');

% Normalized angle handles profile direction internally. Sparsity is zero
% because beta is signed; default 90%% sparsity would discard negative beta.
gm = GradientMaps( ...
    'kernel',cfg.brainspace_kernel, ...
    'approach',cfg.brainspace_approach, ...
    'alignment',cfg.brainspace_alignment, ...
    'n_components',cfg.gradient_components, ...
    'random_state',cfg.gradient_random_state, ...
    'verbose',true);
gm = gm.fit(profile, ...
    'sparsity',cfg.brainspace_sparsity, ...
    'alpha',cfg.brainspace_alpha, ...
    'diffusion_time',cfg.brainspace_diffusion_time);

gradients = single(gm.gradients{1});
eigenvalues = double(gm.lambda{1}(:));
brainspace_method = gm.method;

% BrainSpace gradient signs are arbitrary. Preserve the raw orientation.
for g = 1:size(gradients,2)
    write_map_like(gradients(:,g),wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,sprintf('group_brainspace_gradient-%02d.nii.gz',g)));
end

parameter = ["kernel";"approach";"alignment";"sparsity";"alpha"; ...
    "diffusion_time";"n_components";"n_voxels";"n_features"];
value = [string(cfg.brainspace_kernel);string(cfg.brainspace_approach); ...
    string(cfg.brainspace_alignment);string(cfg.brainspace_sparsity); ...
    string(cfg.brainspace_alpha);string(cfg.brainspace_diffusion_time); ...
    string(cfg.gradient_components);string(size(profile,1)); ...
    string(size(profile,2))];
writetable(table(parameter,value), ...
    fullfile(out_dir,'brainspace_parameter_record.csv'));
save(fullfile(out_dir,'group_brainspace_gradient_results.mat'), ...
    'gradients','eigenvalues','wm_idx', ...
    'brainspace_method','-v7.3');

fig = figure('Color','w','Position',[50 50 1400 850],'Visible','off');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile;
n_show = min(20,numel(eigenvalues));
plot(1:n_show,eigenvalues(1:n_show),'o-','LineWidth',1.5);
xlabel('Component'); ylabel('Eigenvalue');
title('BrainSpace diffusion eigenspectrum'); grid on;
nexttile;
relative = eigenvalues(1:n_show)/sum(max(eigenvalues(1:n_show),0));
bar(1:n_show,100*relative);
xlabel('Component'); ylabel('Relative eigenvalue (%)');
title('Relative eigenspectrum');
nexttile;
histogram(gradients(:,1),100);
xlabel('Gradient 1 score'); ylabel('WM voxels');
title('Gradient 1 distribution');

z_slices = [38 46 54];
mask = niftiread(cfg.wm_mask_final)>0;
vol = zeros(size(mask),'single');
vol(wm_idx) = gradients(:,1);
lim = max(abs(prctile(double(gradients(:,1)),[1 99])));
for j = 1:3
    nexttile;
    img = rot90(vol(:,:,z_slices(j)));
    m = rot90(mask(:,:,z_slices(j)));
    img(~m) = NaN;
    imagesc(img,[-lim lim]); axis image off;
    colormap(gca,turbo(256)); colorbar;
    title(sprintf('Gradient 1, axial z=%d',z_slices(j)));
end
exportgraphics(fig,fullfile(out_dir,'brainspace_gradient_summary.png'), ...
    'Resolution',200);
close(fig);
fprintf('BrainSpace gradients saved to: %s\n',out_dir);
