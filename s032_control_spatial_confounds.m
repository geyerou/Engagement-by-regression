%% s032_control_spatial_confounds
% Residualize key local maps against anatomical/signal/spatial confounds.
% The goal is not to "explain away" anatomy, but to test whether the local
% effects are only trivial GM proximity / WM-boundary / signal-quality maps.
clear; clc; [cfg,ext]=round1_local_extension_config();
out_dir=ext.result_dirs{12};

C=load(fullfile(ext.result_dirs{11},'spatial_confound_maps.mat'), ...
    'confounds','wm_idx');
L=load(fullfile(ext.result_dirs{1},'local_r2_fc_group_data.mat'), ...
    'R2_residual','zR2','zFC_rms','group_mean_residual','wm_idx');
G=load(fullfile(ext.result_dirs{4},'gradient_subject_inference.mat'), ...
    'difference','mean_difference');
if ~isequal(C.wm_idx,L.wm_idx), error('WM mismatch for local R2 data.'); end

X=build_confound_design(C.confounds);
n=size(L.R2_residual,1); V=numel(C.wm_idx);

R2_resid_controlled=zeros(n,V,'single');
G1diff_controlled=zeros(n,V,'single');
for s=1:n
    R2_resid_controlled(s,:)=single(residualize_vector( ...
        double(L.R2_residual(s,:))',X));
    G1diff_controlled(s,:)=single(residualize_vector( ...
        double(G.difference(s,:))',X));
end
group_R2_residual_controlled=mean(R2_resid_controlled,1)';
group_G1diff_controlled=mean(G1diff_controlled,1)';

% Also create a fully controlled version of the raw zR2-zFC relation:
% zR2 ~ zFC + confounds.  This asks whether local prediction exceeds local FC
% after the same spatial nuisance controls.
R2_minus_FC_confounds=zeros(n,V,'single');
for s=1:n
    Xfull=[ones(V,1), zscore(double(L.zFC_rms(s,:))'), X(:,2:end)];
    R2_minus_FC_confounds(s,:)=single(residualize_vector( ...
        zscore(double(L.zR2(s,:))'),Xfull));
end
group_R2_minus_FC_confounds=mean(R2_minus_FC_confounds,1)';

[r2_t,r2_p,r2_q,r2_pfwer]=signflip_maxT(double(R2_resid_controlled), ...
    ext.n_signflip,20260623);
[g1_t,g1_p,g1_q,g1_pfwer]=signflip_maxT(double(G1diff_controlled), ...
    ext.n_signflip,20260624);
[rf_t,rf_p,rf_q,rf_pfwer]=signflip_maxT(double(R2_minus_FC_confounds), ...
    ext.n_signflip,20260625);

write_map_like(group_R2_residual_controlled,C.wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_R2_residual_after_spatial_confounds.nii.gz'));
write_map_like(group_G1diff_controlled,C.wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_beta_minus_FC_G1_after_spatial_confounds.nii.gz'));
write_map_like(group_R2_minus_FC_confounds,C.wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'group_zR2_after_zFC_and_spatial_confounds.nii.gz'));
write_map_like(single(r2_pfwer),C.wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'R2_residual_controlled_maxT_FWER_p.nii.gz'));
write_map_like(single(g1_pfwer),C.wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'beta_minus_FC_G1_controlled_maxT_FWER_p.nii.gz'));
write_map_like(single(rf_pfwer),C.wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'zR2_after_zFC_spatial_confounds_maxT_FWER_p.nii.gz'));

map_r_before=corr(double(L.group_mean_residual),double(G.mean_difference));
map_r_after=corr(double(group_R2_residual_controlled), ...
    double(group_G1diff_controlled));
joint_before=nnz((abs(zscore(double(L.group_mean_residual)))>2) & ...
    (abs(zscore(double(G.mean_difference)))>2));
joint_after=nnz((r2_pfwer(:)<.05) & (g1_pfwer(:)<.05));
if joint_after>0
    joint_mask=(r2_pfwer(:)<.05) & (g1_pfwer(:)<.05);
    joint_after_sign_agreement=mean(sign(group_R2_residual_controlled(joint_mask)) == ...
        sign(group_G1diff_controlled(joint_mask)));
else
    joint_after_sign_agreement=NaN;
end

metric=["map_r_R2resid_G1diff_before_control"; ...
    "map_r_R2resid_G1diff_after_control"; ...
    "R2_residual_controlled_FWER_voxels"; ...
    "G1diff_controlled_FWER_voxels"; ...
    "zR2_after_zFC_spatial_confounds_FWER_voxels"; ...
    "joint_controlled_FWER_voxels"; ...
    "joint_controlled_sign_agreement"; ...
    "rough_joint_voxels_before_control_z2"];
value=[map_r_before;map_r_after;sum(r2_pfwer<.05);sum(g1_pfwer<.05); ...
    sum(rf_pfwer<.05);joint_after;joint_after_sign_agreement;joint_before];
writetable(table(metric,value),fullfile(out_dir, ...
    'spatial_confound_control_summary.csv'));

save(fullfile(out_dir,'spatial_confound_control.mat'), ...
    'R2_resid_controlled','G1diff_controlled','R2_minus_FC_confounds', ...
    'group_R2_residual_controlled','group_G1diff_controlled', ...
    'group_R2_minus_FC_confounds','r2_t','r2_p','r2_q','r2_pfwer', ...
    'g1_t','g1_p','g1_q','g1_pfwer','rf_t','rf_p','rf_q','rf_pfwer', ...
    'X','C','map_r_before','map_r_after','-v7.3');
fprintf('s032 complete: controlled map r %.3f -> %.3f\n', ...
    map_r_before,map_r_after);

function X=build_confound_design(T)
dgm=double(T.dist_GM_mm); dwm=double(T.dist_WM_boundary_mm);
prev=double(T.WM_prevalence); sd=double(T.temporal_SD);
x=double(T.x_mm); y=double(T.y_mm); z=double(T.z_mm);
cols=[log1p(dgm),log1p(dwm),prev,log1p(sd),x,y,z,x.^2,y.^2,z.^2, ...
    x.*y,x.*z,y.*z];
if any(isfinite(double(T.dist_CSF_mm)))
    cols=[cols,log1p(double(T.dist_CSF_mm))]; %#ok<AGROW>
end
bad=all(~isfinite(cols),1) | std(cols,0,1,'omitnan')<eps;
cols(:,bad)=[];
for c=1:size(cols,2)
    v=cols(:,c); v(~isfinite(v))=median(v,'omitnan');
    cols(:,c)=zscore(v);
end
X=[ones(height(T),1),cols];
end

function r=residualize_vector(y,X)
y=y(:); y(~isfinite(y))=median(y,'omitnan');
b=X\y; r=y-X*b; r=zscore(r);
end

function [t_obs,p,q,p_fwer]=signflip_maxT(X,nperm,seed)
n=size(X,1); V=size(X,2);
t_obs=mean(X,1)./max(std(X,0,1)/sqrt(n),eps);
p=2*tcdf(-abs(t_obs),n-1); q=fdr_bh(p);
rng(seed,'twister'); max_abs_t=zeros(nperm,1);
for k=1:nperm
    signs=2*(rand(n,1)>.5)-1;
    P=X.*signs;
    tp=mean(P,1)./max(std(P,0,1)/sqrt(n),eps);
    max_abs_t(k)=max(abs(tp));
end
p_fwer=zeros(1,V);
for first=1:5000:V
    ii=first:min(first+4999,V);
    p_fwer(ii)=(1+sum(max_abs_t>=abs(t_obs(ii)),1))/(1+nperm);
end
end
