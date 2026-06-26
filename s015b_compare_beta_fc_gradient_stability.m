%% s015b_compare_beta_fc_gradient_stability
% Compare split-half BrainSpace gradient stability for group beta and FC
% fingerprints. Run s015a first.
%
% For each split scheme and modality, BrainSpace jointly fits the two split
% profiles and performs Procrustes alignment. This is memory-intensive:
% run this script alone.

clear; clc;
cfg = load_project_config();
if ~exist(cfg.brainspace_path,'dir')
    error('BrainSpace path not found: %s',cfg.brainspace_path);
end
addpath(genpath(cfg.brainspace_path));

profile_dir = fullfile(cfg.result_dirs{16},'s015a_beta_fc_stability');
out_dir = fullfile(cfg.result_dirs{16},'s015b_gradient_stability');
ensure_dir(out_dir);
split_names = string(cfg.stability_split_names);
modalities = ["beta","FC"];
summary_rows = cell(numel(split_names)*numel(modalities), ...
    5+cfg.gradient_components);
row_counter = 0;

for q = 1:numel(split_names)
    split_name = split_names(q);
    profile_file = fullfile(profile_dir,sprintf( ...
        'group_%s_split_profiles.mat',split_name));
    if ~exist(profile_file,'file')
        error('Missing s015a split profile: %s',profile_file);
    end
    D = load(profile_file);

    for m = 1:numel(modalities)
        modality = modalities(m);
        if modality=="beta"
            profile_a = double(D.group_beta_split_a);
            profile_b = double(D.group_beta_split_b);
        else
            profile_a = double(D.group_fc_split_a);
            profile_b = double(D.group_fc_split_b);
        end
        if any(~isfinite(profile_a(:))) || any(~isfinite(profile_b(:)))
            error('%s %s profiles contain NaN/Inf.',split_name,modality);
        end

        fprintf('BrainSpace stability: %s, %s\n',split_name,modality);
        gm = GradientMaps( ...
            'kernel',cfg.brainspace_kernel, ...
            'approach',cfg.brainspace_approach, ...
            'alignment','pa', ...
            'n_components',cfg.gradient_components, ...
            'random_state',cfg.gradient_random_state, ...
            'verbose',true);
        gm = gm.fit({profile_a,profile_b}, ...
            'sparsity',cfg.brainspace_sparsity, ...
            'alpha',cfg.brainspace_alpha, ...
            'diffusion_time',cfg.brainspace_diffusion_time);

        raw_a = single(gm.gradients{1});
        raw_b = single(gm.gradients{2});
        aligned_a = single(gm.aligned{1});
        aligned_b = single(gm.aligned{2});
        eigenvalues_a = double(gm.lambda{1}(:));
        eigenvalues_b = double(gm.lambda{2}(:));

        component_r = zeros(1,cfg.gradient_components);
        for g = 1:cfg.gradient_components
            component_r(g) = corr(double(aligned_a(:,g)), ...
                double(aligned_b(:,g)));
            write_map_like(aligned_a(:,g),D.wm_idx,cfg.wm_mask_final, ...
                fullfile(out_dir,sprintf( ...
                '%s_%s_splitA_aligned_gradient-%02d.nii.gz', ...
                split_name,modality,g)));
            write_map_like(aligned_b(:,g),D.wm_idx,cfg.wm_mask_final, ...
                fullfile(out_dir,sprintf( ...
                '%s_%s_splitB_aligned_gradient-%02d.nii.gz', ...
                split_name,modality,g)));
        end

        row_counter = row_counter+1;
        summary_rows(row_counter,:) = [{split_name,modality, ...
            mean(component_r(1:min(3,end))),component_r(1), ...
            component_r(2)},num2cell(component_r)];

        save(fullfile(out_dir,sprintf( ...
            '%s_%s_gradient_stability.mat',split_name,modality)), ...
            'raw_a','raw_b','aligned_a','aligned_b','component_r', ...
            'eigenvalues_a','eigenvalues_b','-v7.3');
        clear profile_a profile_b raw_a raw_b aligned_a aligned_b gm
    end
end

variable_names = [{'split_name','modality','mean_r_first3', ...
    'gradient1_r','gradient2_r'}, ...
    arrayfun(@(g)sprintf('gradient%02d_r',g), ...
    1:cfg.gradient_components,'UniformOutput',false)];
summary = cell2table(summary_rows,'VariableNames',variable_names);
writetable(summary,fullfile(out_dir,'beta_vs_FC_gradient_stability.csv'));

fig = figure('Color','w','Position',[50 50 1100 500],'Visible','off');
for q = 1:numel(split_names)
    subplot(1,numel(split_names),q);
    rows = summary.split_name==split_names(q);
    vals = [summary.gradient01_r(rows),summary.gradient02_r(rows), ...
        summary.gradient03_r(rows)];
    bar(vals');
    ylim([-1 1]); grid on;
    xlabel('Gradient'); ylabel('Split-half correlation');
    title(sprintf('%s split',split_names(q)));
    legend(cellstr(summary.modality(rows)),'Location','best');
end
exportgraphics(fig,fullfile(out_dir, ...
    'beta_vs_FC_gradient_stability.png'),'Resolution',200);
close(fig);

fprintf('Gradient stability results saved to: %s\n',out_dir);
