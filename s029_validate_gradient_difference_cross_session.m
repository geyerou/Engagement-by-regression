%% s029_validate_gradient_difference_cross_session
% Cross-session validation: projections learned from one group session are
% applied only to subject fingerprints from the opposite session.
clear; clc; [cfg,ext]=round1_local_extension_config();
profile_dir=fullfile(cfg.result_dirs{16},'s015a_beta_fc_stability');
gradient_dir=fullfile(cfg.result_dirs{16},'s015b_gradient_stability');
D=load(fullfile(profile_dir,'group_session_split_profiles.mat'));
ba=double(niftiread(fullfile(gradient_dir, ...
    'session_beta_splitA_aligned_gradient-01.nii.gz')));
bb=double(niftiread(fullfile(gradient_dir, ...
    'session_beta_splitB_aligned_gradient-01.nii.gz')));
fa=double(niftiread(fullfile(gradient_dir, ...
    'session_FC_splitA_aligned_gradient-01.nii.gz')));
fb=double(niftiread(fullfile(gradient_dir, ...
    'session_FC_splitB_aligned_gradient-01.nii.gz')));
targets={zscore(ba(D.wm_idx)),zscore(bb(D.wm_idx)); ...
    zscore(fa(D.wm_idx)),zscore(fb(D.wm_idx))};
if corr(targets{1,1},targets{1,2})<0,targets{1,2}=-targets{1,2};end
if corr(targets{2,1},targets{2,2})<0,targets{2,2}=-targets{2,2};end
[bp{1},bmu{1},bsd{1}]=fit_projection(double(D.group_beta_split_a), ...
    targets{1,1},ext.gradient_projection_lambda);
[bp{2},bmu{2},bsd{2}]=fit_projection(double(D.group_beta_split_b), ...
    targets{1,2},ext.gradient_projection_lambda);
[fp{1},fmu{1},fsd{1}]=fit_projection(double(D.group_fc_split_a), ...
    targets{2,1},ext.gradient_projection_lambda);
[fp{2},fmu{2},fsd{2}]=fit_projection(double(D.group_fc_split_b), ...
    targets{2,2},ext.gradient_projection_lambda);
L=load(fullfile(cfg.result_dirs{5},'shared_lambda_refit', ...
    'shared_lambda_definition.mat'),'shared_lambda');
validation_lambda=L.shared_lambda;

n=numel(cfg.subject_list); V=numel(D.wm_idx);
difference=zeros(n,V,2,'single'); direction_r=zeros(n,1);
for s=1:n
    sid=cfg.subject_list{s};
    fprintf('s029 subject %s (%d/%d)\n',sid,s,n);
    [X_runs,Y_runs,meta]=load_subject_runs(fullfile(cfg.result_dirs{4}, ...
        sprintf('sub-%s_model_input.mat',sid)));
    if ~isequal(meta.wm_voxel_indices,D.wm_idx),error('WM mismatch.');end
    profiles=cell(2,2);
    for ses=1:2
        rr=(ses-1)*2+(1:2); X=vertcat(X_runs{rr}); Y=vertcat(Y_runs{rr});
        B=fit_final_beta_blocked(X,Y,validation_lambda,cfg.voxel_block_size);
        profiles{1,ses}=double(B'); profiles{2,ses}=double( ...
            compute_fc_profile(X,Y,cfg.voxel_block_size));
    end
    % Direction 1: train REST1 projection, apply REST2 fingerprints.
    zb=apply_projection(profiles{1,2},bmu{1},bsd{1},bp{1});
    zf=apply_projection(profiles{2,2},fmu{1},fsd{1},fp{1});
    difference(s,:,1)=single(zscore(zb)-zscore(zf));
    % Direction 2: train REST2 projection, apply REST1 fingerprints.
    zb=apply_projection(profiles{1,1},bmu{2},bsd{2},bp{2});
    zf=apply_projection(profiles{2,1},fmu{2},fsd{2},fp{2});
    difference(s,:,2)=single(zscore(zb)-zscore(zf));
    direction_r(s)=corr(double(difference(s,:,1))', ...
        double(difference(s,:,2))');
    clear X_runs Y_runs profiles
end

mean_difference=squeeze(mean(difference,1));
group_map_r=corr(double(mean_difference(:,1)),double(mean_difference(:,2)));
crossfit_mean=mean(mean_difference,2);
X=mean(double(difference),3);
[t_obs,p_fwer]=max_t_signflip(X,ext.n_signflip,20260629);
write_map_like(crossfit_mean,D.wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{9}, ...
    'cross_session_mean_beta_minus_FC_G1.nii.gz'));
write_map_like(single(p_fwer),D.wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{9}, ...
    'cross_session_beta_minus_FC_G1_maxT_FWER_p.nii.gz'));
subject_id=string(cfg.subject_list(:));
writetable(table(subject_id,direction_r,'VariableNames', ...
    {'subject_id','REST1_to_REST2_vs_REST2_to_REST1_spatial_r'}), ...
    fullfile(ext.result_dirs{9}, ...
    'cross_session_gradient_subject_summary.csv'));
metric=["mean_subject_cross_direction_r";"median_subject_cross_direction_r"; ...
    "group_cross_direction_map_r";"cross_session_maxT_voxels"];
value=[mean(direction_r);median(direction_r);group_map_r;sum(p_fwer<.05)];
writetable(table(metric,value),fullfile(ext.result_dirs{9}, ...
    'cross_session_gradient_group_summary.csv'));
save(fullfile(ext.result_dirs{9},'cross_session_gradient_validation.mat'), ...
    'difference','direction_r','mean_difference','crossfit_mean', ...
    'group_map_r','t_obs','p_fwer','validation_lambda','-v7.3');
fprintf('s029 complete: mean subject r=%.3f; group map r=%.3f; maxT=%d.\n', ...
    mean(direction_r),group_map_r,sum(p_fwer<.05));

function [b,mu,sd]=fit_projection(X,y,lambda)
mu=mean(X,1); sd=std(X,0,1); sd(sd<eps)=1; Z=(X-mu)./sd;
b=(Z'*Z+size(Z,1)*lambda*eye(size(Z,2)))\(Z'*y);
end
function score=apply_projection(X,mu,sd,b)
score=((X-mu)./sd)*b;
end
function [t_obs,p_fwer]=max_t_signflip(X,n_perm,seed)
[n,V]=size(X); t_obs=mean(X,1)./max(std(X,0,1)/sqrt(n),eps);
rng(seed,'twister'); mx=zeros(n_perm,1);
for p=1:n_perm
    sg=2*(rand(n,1)>.5)-1; P=X.*sg;
    tp=mean(P,1)./max(std(P,0,1)/sqrt(n),eps); mx(p)=max(abs(tp));
end
p_fwer=zeros(1,V);
for first=1:5000:V
    ii=first:min(first+4999,V);
    p_fwer(ii)=(1+sum(mx>=abs(t_obs(ii)),1))/(1+n_perm);
end
end
