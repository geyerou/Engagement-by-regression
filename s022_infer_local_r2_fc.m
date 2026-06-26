%% s022_infer_local_r2_fc
% Voxelwise max-|t| FWER and signed cluster-mass sign-flip inference.
% Positive and negative clusters are formed separately. The permutation
% maximum is taken across both signs, controlling cluster FWER jointly.
clear; clc; [cfg,ext]=round1_local_extension_config();
D=load(fullfile(ext.result_dirs{1},'local_r2_fc_group_data.mat'), ...
    'R2_residual','wm_idx');
X=double(D.R2_residual); [n,V]=size(X);
t_obs=mean(X,1)./max(std(X,0,1)/sqrt(n),eps);
tcrit=tinv(1-ext.cluster_forming_p/2,n-1);
wm=niftiread(cfg.wm_mask_final)>0;
obs_vol=zeros(size(wm)); obs_vol(D.wm_idx)=t_obs;
[CC_pos,CC_neg,obs_mass_pos,obs_mass_neg]=signed_clusters( ...
    obs_vol,wm,tcrit);

rng(20260620,'twister');
max_abs_t=zeros(ext.n_signflip,1);
max_cluster_mass=zeros(ext.n_signflip,1);
for p=1:ext.n_signflip
    signs=2*(rand(n,1)>.5)-1; P=X.*signs;
    tp=mean(P,1)./max(std(P,0,1)/sqrt(n),eps);
    max_abs_t(p)=max(abs(tp));
    vol=zeros(size(wm)); vol(D.wm_idx)=tp;
    [~,~,mass_pos,mass_neg]=signed_clusters(vol,wm,tcrit);
    max_cluster_mass(p)=max([0;mass_pos(:);mass_neg(:)]);
end
p_fwer=zeros(1,V);
for first=1:5000:V
    ii=first:min(first+4999,V);
    p_fwer(ii)=(1+sum(max_abs_t>=abs(t_obs(ii)),1))/ ...
        (1+ext.n_signflip);
end
cluster_map=zeros(size(wm),'single');
cluster_id=zeros(size(wm),'single');
pixel_lists=[CC_pos.PixelIdxList,CC_neg.PixelIdxList];
direction=[repmat("positive",CC_pos.NumObjects,1); ...
    repmat("negative",CC_neg.NumObjects,1)];
obs_mass=[obs_mass_pos(:);obs_mass_neg(:)];
cluster_p=ones(numel(pixel_lists),1);
for c=1:numel(pixel_lists)
    cluster_p(c)=(1+sum(max_cluster_mass>=obs_mass(c)))/(1+ext.n_signflip);
    signed_id=c;
    if direction(c)=="negative",signed_id=-c;end
    cluster_id(pixel_lists{c})=signed_id;
    cluster_map(pixel_lists{c})=cluster_p(c);
end
write_map_like(single(t_obs),D.wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{2},'local_R2_residual_t.nii.gz'));
write_map_like(single(p_fwer),D.wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{2},'local_R2_residual_maxT_FWER_p.nii.gz'));
niftiwrite(cluster_map,fullfile(ext.result_dirs{2}, ...
    'local_R2_residual_cluster_mass_p.nii'), ...
    niftiinfo(cfg.wm_mask_final),'Compressed',true);
niftiwrite(cluster_id,fullfile(ext.result_dirs{2}, ...
    'local_R2_residual_cluster_id.nii'), ...
    niftiinfo(cfg.wm_mask_final),'Compressed',true);
cluster_table=table((1:numel(pixel_lists))',direction, ...
    cellfun(@numel,pixel_lists)',obs_mass,cluster_p,'VariableNames', ...
    {'cluster_id','direction','voxels','signed_t_mass','permutation_p'});
writetable(cluster_table,fullfile(ext.result_dirs{2}, ...
    'local_R2_residual_clusters.csv'));
save(fullfile(ext.result_dirs{2},'local_r2_fc_inference.mat'), ...
    't_obs','p_fwer','cluster_table','max_abs_t','max_cluster_mass', ...
    'tcrit','-v7.3');
fprintf('s022 complete: %d clusters; %d cluster-mass FWER p<.05.\n', ...
    height(cluster_table),sum(cluster_p<.05));

function [CC_pos,CC_neg,mass_pos,mass_neg]=signed_clusters(vol,wm,tcrit)
CC_pos=bwconncomp(wm & vol>=tcrit,26);
CC_neg=bwconncomp(wm & vol<=-tcrit,26);
mass_pos=cellfun(@(ii)sum(vol(ii)),CC_pos.PixelIdxList);
mass_neg=cellfun(@(ii)sum(-vol(ii)),CC_neg.PixelIdxList);
mass_pos=mass_pos(:); mass_neg=mass_neg(:);
end
