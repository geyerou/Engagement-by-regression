%% s024_infer_subject_beta_fc_gradient_difference
% Project each subject's fingerprints onto fixed beta and FC G1 templates.
clear; clc; [cfg,ext]=round1_local_extension_config();
T=load(fullfile(ext.result_dirs{3}, ...
    'gradient_templates_and_projections.mat'));
n=numel(cfg.subject_list); V=numel(T.wm_idx);
beta_score=zeros(n,V,'single'); fc_score=zeros(n,V,'single');
difference=zeros(n,V,'single'); template_r=zeros(n,2);
for s=1:n
    sid=cfg.subject_list{s};
    B=load(fullfile(cfg.result_dirs{7},sprintf( ...
        'sub-%s_beta_profile.mat',sid)),'Beta_profile','wm_voxel_indices');
    F=load(fullfile(cfg.result_dirs{12},sprintf( ...
        'sub-%s_FC_profile.mat',sid)),'FC_profile','wm_voxel_indices');
    if ~isequal(B.wm_voxel_indices,T.wm_idx) || ...
            ~isequal(F.wm_voxel_indices,T.wm_idx)
        error('WM mismatch for %s.',sid);
    end
    zb=((double(B.Beta_profile)-T.beta_mu)./T.beta_sd)*T.beta_projection;
    zf=((double(F.FC_profile)-T.fc_mu)./T.fc_sd)*T.fc_projection;
    zb=zscore(zb); zf=zscore(zf);
    beta_score(s,:)=single(zb); fc_score(s,:)=single(zf);
    difference(s,:)=single(zb-zf);
    template_r(s,:)=[corr(zb,T.beta_template),corr(zf,T.fc_template)];
    save(fullfile(ext.result_dirs{4},sprintf( ...
        'sub-%s_beta_FC_G1_scores.mat',sid)),'sid','zb','zf','-v7.3');
end

X=double(difference); t_obs=mean(X,1)./max(std(X,0,1)/sqrt(n),eps);
p=2*tcdf(-abs(t_obs),n-1); q=fdr_bh(p);
rng(20260620,'twister'); max_abs_t=zeros(ext.n_signflip,1);
for k=1:ext.n_signflip
    signs=2*(rand(n,1)>.5)-1; P=X.*signs;
    tp=mean(P,1)./max(std(P,0,1)/sqrt(n),eps);
    max_abs_t(k)=max(abs(tp));
end
p_fwer=zeros(1,V);
for first=1:5000:V
    ii=first:min(first+4999,V);
    p_fwer(ii)=(1+sum(max_abs_t>=abs(t_obs(ii)),1))/ ...
        (1+ext.n_signflip);
end
mean_difference=mean(difference,1)';
write_map_like(mean_difference,T.wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{4},'subject_mean_beta_minus_FC_G1.nii.gz'));
write_map_like(single(q),T.wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{4},'beta_minus_FC_G1_FDR_q.nii.gz'));
write_map_like(single(p_fwer),T.wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{4},'beta_minus_FC_G1_maxT_FWER_p.nii.gz'));
subject_id=string(cfg.subject_list(:));
subject_summary=table(subject_id,template_r(:,1),template_r(:,2), ...
    mean(abs(difference),2),'VariableNames',{'subject_id', ...
    'beta_template_r','FC_template_r','mean_absolute_G1_difference'});
writetable(subject_summary,fullfile(ext.result_dirs{4}, ...
    'gradient_subject_summary.csv'));
save(fullfile(ext.result_dirs{4},'gradient_subject_inference.mat'), ...
    'beta_score','fc_score','difference','mean_difference','t_obs', ...
    'p','q','p_fwer','max_abs_t','template_r','T','-v7.3');
fprintf('s024 complete: FDR voxels=%d; maxT-FWER voxels=%d.\n', ...
    sum(q<.05),sum(p_fwer<.05));
