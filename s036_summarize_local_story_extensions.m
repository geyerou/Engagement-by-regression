%% s036_summarize_local_story_extensions
% One-stop summary for s031-s035 extensions.
clear; clc; [cfg,ext]=round1_local_extension_config();
out_dir=ext.result_dirs{16};

C=read_metric_table(fullfile(ext.result_dirs{11}, ...
    'spatial_confound_summary.csv'));
K=read_metric_table(fullfile(ext.result_dirs{12}, ...
    'spatial_confound_control_summary.csv'));
N=read_metric_table(fullfile(ext.result_dirs{13}, ...
    'spatial_block_null_summary.csv'));
B=readtable(fullfile(ext.result_dirs{14},'G1_bin_summary.csv'), ...
    'Delimiter',',','TextType','string','VariableNamingRule','preserve');
H=read_metric_table(fullfile(ext.result_dirs{15}, ...
    'hemispheric_symmetry_summary.csv'));
HM=readtable(fullfile(ext.result_dirs{15},'hemispheric_symmetry_map_correlations.csv'), ...
    'Delimiter',',','TextType','string','VariableNamingRule','preserve');

metric=[prepend_metric(C.metric,"confound_"); ...
    prepend_metric(K.metric,"control_"); ...
    prepend_metric(N.metric,"spatial_null_"); ...
    prepend_metric(H.metric,"symmetry_")];
value=[C.value;K.value;N.value;H.value];
writetable(table(metric,value),fullfile(out_dir, ...
    'local_story_extension_summary.csv'));

fig=figure('Color','w','Visible','off','Position',[50 50 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
bar([getv(K,'map_r_R2resid_G1diff_before_control'), ...
    getv(K,'map_r_R2resid_G1diff_after_control')]);
xticklabels({'Before','After'}); ylim([-1 1]); ylabel('Spatial r');
title('R2 residual vs beta-FC G1 after confound control');

nexttile;
bar([getv(K,'R2_residual_controlled_FWER_voxels'), ...
    getv(K,'G1diff_controlled_FWER_voxels'), ...
    getv(K,'joint_controlled_FWER_voxels')]);
xticklabels({'R2 resid','G1 diff','Joint'}); ylabel('FWER voxels');
title('Controlled significant voxels');

nexttile;
bar([getv(N,'controlled_R2resid_G1diff_block_p'), ...
    getv(N,'controlled_same_sign_joint_block_p')]);
xticklabels({'Map r','Same-sign joint'}); ylim([0 1]); ylabel('Block-null p');
title('Spatial block permutation tests');

nexttile;
plot(B.bin,B.mean_R2_residual,'-o'); hold on;
plot(B.bin,B.mean_R2_residual_controlled,'-o');
plot(B.bin,B.mean_zR2_after_zFC_spatial_confounds,'-o');
legend({'R2 residual','controlled R2 residual','zR2|zFC+conf'}, ...
    'Location','best'); xlabel('beta-FC G1 bin');
title('Prediction structure along G1');

nexttile;
plot(B.bin,B.mean_dist_GM_mm,'-o'); hold on;
plot(B.bin,B.mean_WM_prevalence,'-o');
yyaxis right; plot(B.bin,B.mean_zFC_RMS,'-o');
xlabel('beta-FC G1 bin'); title('Contamination/signal diagnostics');
legend({'GM distance','WM prevalence','zFC RMS'},'Location','best');

nexttile;
bar(categorical(HM.map_name),HM.LR_correlation); ylim([-1 1]);
xtickangle(35); ylabel('LR mirror r'); title('Hemispheric symmetry');

exportgraphics(fig,fullfile(out_dir,'local_story_extension_summary.png'), ...
    'Resolution',200);
close(fig);

copyfile(fullfile(ext.result_dirs{14},'G1_bin_summary.csv'), ...
    fullfile(out_dir,'G1_bin_summary.csv'));
copyfile(fullfile(ext.result_dirs{14},'G1_bin_network_contribution.csv'), ...
    fullfile(out_dir,'G1_bin_network_contribution.csv'));
copyfile(fullfile(ext.result_dirs{15},'hemispheric_symmetry_map_correlations.csv'), ...
    fullfile(out_dir,'hemispheric_symmetry_map_correlations.csv'));

fprintf('s036 complete: story summary saved to %s\n',out_dir);

function out=prepend_metric(metric,prefix)
out=string(metric);
for i=1:numel(out)
    out(i)=prefix+out(i);
end
end

function v=getv(T,name)
idx=string(T.metric)==string(name);
if ~any(idx), v=NaN; else, v=T.value(find(idx,1)); end
end

function T=read_metric_table(file)
R=readtable(file,'Delimiter',',','TextType','string', ...
    'VariableNamingRule','preserve');
metric=string(R{:,1});
value=R{:,2};
if ~isnumeric(value)
    value=str2double(string(value));
end
T=table(metric,double(value),'VariableNames',{'metric','value'});
end
