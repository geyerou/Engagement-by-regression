%% s033_spatial_autocorr_block_null
% Spatial-autocorrelation-aware block permutation null for key WM map
% relationships.  This is a volumetric WM alternative to cortical spin tests.
clear; clc; [cfg,ext]=round1_local_extension_config();
out_dir=ext.result_dirs{13};

C=load(fullfile(ext.result_dirs{11},'spatial_confound_maps.mat'), ...
    'confounds','wm_idx');
L=load(fullfile(ext.result_dirs{1},'local_r2_fc_group_data.mat'), ...
    'group_mean_residual');
G=load(fullfile(ext.result_dirs{4},'gradient_subject_inference.mat'), ...
    'mean_difference');
K=load(fullfile(ext.result_dirs{12},'spatial_confound_control.mat'), ...
    'group_R2_residual_controlled','group_G1diff_controlled', ...
    'r2_pfwer','g1_pfwer');

x=double(C.confounds.x_mm); y=double(C.confounds.y_mm); z=double(C.confounds.z_mm);
block_id=make_spatial_blocks(x,y,z,ext.spatial_null_block_mm);
nblock=max(block_id); nperm=ext.spatial_null_n_permutations;

raw_r2=double(L.group_mean_residual(:));
raw_g1=double(G.mean_difference(:));
ctl_r2=double(K.group_R2_residual_controlled(:));
ctl_g1=double(K.group_G1diff_controlled(:));

obs_raw_r=corr(raw_r2,raw_g1);
obs_ctl_r=corr(ctl_r2,ctl_g1);
r2_sig=K.r2_pfwer(:)<.05; g1_sig=K.g1_pfwer(:)<.05;
obs_joint=nnz(r2_sig & g1_sig);
obs_same_sign=nnz(r2_sig & g1_sig & (sign(ctl_r2)==sign(ctl_g1)));

null_raw_r=zeros(nperm,1);
null_ctl_r=zeros(nperm,1);
null_joint=zeros(nperm,1);
null_same_sign=zeros(nperm,1);
rng(20260626,'twister');
for p=1:nperm
    perm=randperm(nblock);
    raw_g1_perm=permute_by_block(raw_g1,block_id,perm);
    ctl_g1_perm=permute_by_block(ctl_g1,block_id,perm);
    g1_sig_perm=permute_by_block(double(g1_sig),block_id,perm)>.5;
    null_raw_r(p)=corr(raw_r2,raw_g1_perm);
    null_ctl_r(p)=corr(ctl_r2,ctl_g1_perm);
    null_joint(p)=nnz(r2_sig & g1_sig_perm);
    null_same_sign(p)=nnz(r2_sig & g1_sig_perm & ...
        (sign(ctl_r2)==sign(ctl_g1_perm)));
end

p_raw_r=(1+sum(abs(null_raw_r)>=abs(obs_raw_r)))/(1+nperm);
p_ctl_r=(1+sum(abs(null_ctl_r)>=abs(obs_ctl_r)))/(1+nperm);
p_joint=(1+sum(null_joint>=obs_joint))/(1+nperm);
p_same_sign=(1+sum(null_same_sign>=obs_same_sign))/(1+nperm);

metric=["raw_R2resid_G1diff_spatial_r"; ...
    "raw_R2resid_G1diff_block_p"; ...
    "controlled_R2resid_G1diff_spatial_r"; ...
    "controlled_R2resid_G1diff_block_p"; ...
    "controlled_joint_FWER_voxels"; ...
    "controlled_joint_block_p"; ...
    "controlled_same_sign_joint_voxels"; ...
    "controlled_same_sign_joint_block_p"; ...
    "spatial_block_size_mm"; ...
    "n_spatial_blocks"; ...
    "n_permutations"];
value=[obs_raw_r;p_raw_r;obs_ctl_r;p_ctl_r;obs_joint;p_joint; ...
    obs_same_sign;p_same_sign;ext.spatial_null_block_mm;nblock;nperm];
writetable(table(metric,value),fullfile(out_dir, ...
    'spatial_block_null_summary.csv'));

fig=figure('Color','w','Visible','off','Position',[50 50 1200 430]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; histogram(null_raw_r,40); hold on; xline(obs_raw_r,'r','LineWidth',2);
xlabel('Block-null r'); ylabel('Permutations');
title(sprintf('Raw map r p=%.4f',p_raw_r));
nexttile; histogram(null_ctl_r,40); hold on; xline(obs_ctl_r,'r','LineWidth',2);
xlabel('Block-null r'); title(sprintf('Controlled map r p=%.4f',p_ctl_r));
nexttile; histogram(null_same_sign,40); hold on; xline(obs_same_sign,'r','LineWidth',2);
xlabel('Same-sign joint voxels'); title(sprintf('Joint p=%.4f',p_same_sign));
exportgraphics(fig,fullfile(out_dir,'spatial_block_null_summary.png'), ...
    'Resolution',200);
close(fig);

save(fullfile(out_dir,'spatial_block_null.mat'),'block_id','null_raw_r', ...
    'null_ctl_r','null_joint','null_same_sign','obs_raw_r','obs_ctl_r', ...
    'obs_joint','obs_same_sign','p_raw_r','p_ctl_r','p_joint', ...
    'p_same_sign','-v7.3');
fprintf('s033 complete: controlled r %.3f, block p %.4f\n', ...
    obs_ctl_r,p_ctl_r);

function block_id=make_spatial_blocks(x,y,z,block_mm)
bx=floor((x-min(x))/block_mm)+1;
by=floor((y-min(y))/block_mm)+1;
bz=floor((z-min(z))/block_mm)+1;
[~,~,block_id]=unique([bx(:),by(:),bz(:)],'rows');
end

function yp=permute_by_block(y,block_id,perm)
yp=zeros(size(y));
for b=1:max(block_id)
    dst=find(block_id==b);
    src=find(block_id==perm(b));
    if isempty(src)
        yp(dst)=nan;
    elseif numel(src)==numel(dst)
        yp(dst)=y(src);
    else
        pick=src(randi(numel(src),numel(dst),1));
        yp(dst)=y(pick);
    end
end
end
