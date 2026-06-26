%% s035_hemispheric_symmetry_analysis
% Test left-right symmetry of the main WM local maps using mirrored voxel
% coordinates in template space.
clear; clc; [cfg,ext]=round1_local_extension_config();
out_dir=ext.result_dirs{15};

C=load(fullfile(ext.result_dirs{11},'spatial_confound_maps.mat'), ...
    'confounds','wm_idx');
L=load(fullfile(ext.result_dirs{1},'local_r2_fc_group_data.mat'), ...
    'group_mean_zR2','group_mean_zFC_rms','group_mean_residual');
G=load(fullfile(ext.result_dirs{4},'gradient_subject_inference.mat'), ...
    'mean_difference');
K=load(fullfile(ext.result_dirs{12},'spatial_confound_control.mat'), ...
    'group_R2_residual_controlled','group_G1diff_controlled', ...
    'group_R2_minus_FC_confounds');
if isfile(fullfile(ext.result_dirs{14},'G1_binned_cortical_contribution.mat'))
    B=load(fullfile(ext.result_dirs{14},'G1_binned_cortical_contribution.mat'), ...
        'bin_id');
    bin_id=B.bin_id;
else
    bin_id=[];
end

x=double(C.confounds.x_mm); y=double(C.confounds.y_mm); z=double(C.confounds.z_mm);
[left_idx,right_idx,mirror_distance]=find_lr_pairs(x,y,z, ...
    ext.symmetry_pair_tolerance_mm);

map_names=["zR2";"zFC_RMS";"R2_residual";"beta_minus_FC_G1"; ...
    "R2_residual_controlled";"beta_minus_FC_G1_controlled"; ...
    "zR2_after_zFC_spatial_confounds"];
maps=[double(L.group_mean_zR2(:)),double(L.group_mean_zFC_rms(:)), ...
    double(L.group_mean_residual(:)),double(G.mean_difference(:)), ...
    double(K.group_R2_residual_controlled(:)), ...
    double(K.group_G1diff_controlled(:)), ...
    double(K.group_R2_minus_FC_confounds(:))];

rows=cell(numel(map_names),6);
for m=1:numel(map_names)
    lv=maps(left_idx,m); rv=maps(right_idx,m);
    rows(m,:)={map_names(m),numel(left_idx),corr(lv,rv), ...
        mean(abs(lv-rv),'omitnan'),mean(lv-rv,'omitnan'), ...
        corr(abs(lv),abs(rv))};
end
S=cell2table(rows,'VariableNames',{'map_name','n_pairs','LR_correlation', ...
    'mean_absolute_LminusR','mean_signed_LminusR','LR_abs_correlation'});

if ~isempty(bin_id)
    bin_agreement=mean(bin_id(left_idx)==bin_id(right_idx));
    mean_bin_absdiff=mean(abs(bin_id(left_idx)-bin_id(right_idx)));
else
    bin_agreement=NaN; mean_bin_absdiff=NaN;
end
metric=["n_LR_pairs";"median_pair_distance_mm";"G1_bin_exact_agreement"; ...
    "G1_bin_mean_absolute_difference"];
value=[numel(left_idx);median(mirror_distance,'omitnan'); ...
    bin_agreement;mean_bin_absdiff];
writetable(table(metric,value),fullfile(out_dir, ...
    'hemispheric_symmetry_summary.csv'));
writetable(S,fullfile(out_dir,'hemispheric_symmetry_map_correlations.csv'));

pair_table=table(C.wm_idx(left_idx),C.wm_idx(right_idx), ...
    single(mirror_distance),single(x(left_idx)),single(y(left_idx)), ...
    single(z(left_idx)),single(x(right_idx)),single(y(right_idx)), ...
    single(z(right_idx)),'VariableNames',{'left_wm_index','right_wm_index', ...
    'mirror_distance_mm','left_x_mm','left_y_mm','left_z_mm', ...
    'right_x_mm','right_y_mm','right_z_mm'});
writetable(pair_table,fullfile(out_dir,'hemispheric_mirror_pairs.csv'));

fig=figure('Color','w','Visible','off','Position',[50 50 1250 520]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;
bar(categorical(S.map_name),S.LR_correlation); ylim([-1 1]);
xtickangle(35); ylabel('Left-right r'); title('Map symmetry');
nexttile;
scatter(maps(left_idx,4),maps(right_idx,4),4,'filled','MarkerFaceAlpha',.2);
xlabel('Left beta-FC G1'); ylabel('Right mirror beta-FC G1');
title(sprintf('G1 difference LR r=%.3f',S.LR_correlation(S.map_name=="beta_minus_FC_G1")));
nexttile;
scatter(maps(left_idx,5),maps(right_idx,5),4,'filled','MarkerFaceAlpha',.2);
xlabel('Left controlled R2 residual'); ylabel('Right mirror');
title(sprintf('Controlled R2 residual LR r=%.3f', ...
    S.LR_correlation(S.map_name=="R2_residual_controlled")));
exportgraphics(fig,fullfile(out_dir,'hemispheric_symmetry_summary.png'), ...
    'Resolution',200);
close(fig);

save(fullfile(out_dir,'hemispheric_symmetry_analysis.mat'), ...
    'S','pair_table','left_idx','right_idx','mirror_distance','bin_id','-v7.3');
fprintf('s035 complete: %d LR pairs, G1 LR r %.3f\n',numel(left_idx), ...
    S.LR_correlation(S.map_name=="beta_minus_FC_G1"));

function [left_idx,right_idx,dist]=find_lr_pairs(x,y,z,tol)
scale=10; % one decimal millimeter key
keys=strings(numel(x),1);
for i=1:numel(x)
    keys(i)=coord_key(x(i),y(i),z(i),scale);
end
M=containers.Map(cellstr(keys),num2cell(1:numel(keys)));
left=find(x<0);
left_idx=[]; right_idx=[]; dist=[];
for ii=1:numel(left)
    i=left(ii);
    key=coord_key(-x(i),y(i),z(i),scale);
    if isKey(M,key)
        j=M(key);
        d=sqrt((x(j)+x(i))^2+(y(j)-y(i))^2+(z(j)-z(i))^2);
        if d<=tol
            left_idx(end+1,1)=i; %#ok<AGROW>
            right_idx(end+1,1)=j; %#ok<AGROW>
            dist(end+1,1)=d; %#ok<AGROW>
        end
    end
end
end

function key=coord_key(x,y,z,scale)
key=sprintf('%d_%d_%d',round(x*scale),round(y*scale),round(z*scale));
end
