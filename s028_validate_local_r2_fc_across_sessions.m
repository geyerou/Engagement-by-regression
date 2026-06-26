%% s028_validate_local_r2_fc_across_sessions
% Independent REST1 and REST2 estimates of local R2 beyond local FC.
% Within each session, each phase-encoding run predicts the other run.
clear; clc; [cfg,ext]=round1_local_extension_config();
wm=niftiread(cfg.wm_mask_final)>0; wm_idx=find(wm);
n=numel(cfg.subject_list); V=numel(wm_idx);
residual=zeros(n,V,2,'single');
zR2=zeros(n,V,2,'single'); zFC=zeros(n,V,2,'single');
session_spatial_r=zeros(n,1);
L=load(fullfile(cfg.result_dirs{5},'shared_lambda_refit', ...
    'shared_lambda_definition.mat'),'shared_lambda');
validation_lambda=L.shared_lambda;

for s=1:n
    sid=cfg.subject_list{s};
    fprintf('s028 subject %s (%d/%d)\n',sid,s,n);
    [X_runs,Y_runs,meta]=load_subject_runs(fullfile(cfg.result_dirs{4}, ...
        sprintf('sub-%s_model_input.mat',sid)));
    if ~isequal(meta.wm_voxel_indices,wm_idx)
        error('WM mismatch for %s.',sid);
    end
    for ses=1:2
        rr=(ses-1)*2+(1:2);
        sse=zeros(1,V); sst=zeros(1,V);
        for fold=1:2
            te=rr(fold); tr=rr(3-fold);
            pred=ridge_fit_predict_blocked(X_runs{tr},Y_runs{tr}, ...
                X_runs{te},validation_lambda,cfg.voxel_block_size,false);
            y=double(Y_runs{te}); p=double(pred);
            sse=sse+sum((y-p).^2,1);
            sst=sst+sum((y-mean(y,1)).^2,1);
        end
        r2=1-sse./max(sst,eps);
        X=vertcat(X_runs{rr}); Y=vertcat(Y_runs{rr});
        fc=compute_fc_profile(X,Y,cfg.voxel_block_size);
        rmsfc=sqrt(mean(double(fc).^2,2));
        zr=zscore(r2(:)); zf=zscore(rmsfc);
        e=zr-[ones(V,1),zf]*([ones(V,1),zf]\zr);
        zR2(s,:,ses)=single(zr); zFC(s,:,ses)=single(zf);
        residual(s,:,ses)=single(e);
    end
    session_spatial_r(s)=corr(double(residual(s,:,1))', ...
        double(residual(s,:,2))');
    clear X_runs Y_runs
end

mean_residual=squeeze(mean(residual,1));
group_map_r=corr(double(mean_residual(:,1)),double(mean_residual(:,2)));
[t_obs,p_fwer]=max_t_signflip(double(residual),ext.n_signflip,20260628);
conjunction=(p_fwer(:,:,1)<.05 & p_fwer(:,:,2)<.05 & ...
    sign(t_obs(:,:,1))==sign(t_obs(:,:,2)));

for ses=1:2
    write_map_like(mean_residual(:,ses),wm_idx,cfg.wm_mask_final, ...
        fullfile(ext.result_dirs{8},sprintf( ...
        'REST%d_mean_local_R2_residual.nii.gz',ses)));
    write_map_like(single(p_fwer(:,:,ses)),wm_idx,cfg.wm_mask_final, ...
        fullfile(ext.result_dirs{8},sprintf( ...
        'REST%d_local_R2_residual_maxT_FWER_p.nii.gz',ses)));
end
write_map_like(single(conjunction(:)),wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{8}, ...
    'REST1_REST2_local_R2_residual_conjunction.nii.gz'));
subject_id=string(cfg.subject_list(:));
T=table(subject_id,session_spatial_r,'VariableNames', ...
    {'subject_id','REST1_REST2_residual_spatial_r'});
writetable(T,fullfile(ext.result_dirs{8}, ...
    'session_local_R2_residual_subject_summary.csv'));
metric=["mean_subject_session_r";"median_subject_session_r"; ...
    "group_mean_map_r";"same_sign_maxT_conjunction_voxels"];
value=[mean(session_spatial_r);median(session_spatial_r); ...
    group_map_r;nnz(conjunction)];
writetable(table(metric,value),fullfile(ext.result_dirs{8}, ...
    'session_local_R2_residual_group_summary.csv'));
save(fullfile(ext.result_dirs{8},'session_local_R2_fc_validation.mat'), ...
    'residual','zR2','zFC','session_spatial_r','mean_residual', ...
    'group_map_r','t_obs','p_fwer','conjunction','wm_idx', ...
    'validation_lambda','-v7.3');
fprintf('s028 complete: mean subject r=%.3f; group map r=%.3f; conjunction=%d.\n', ...
    mean(session_spatial_r),group_map_r,nnz(conjunction));

function [t_obs,p_fwer]=max_t_signflip(X,n_perm,seed)
[n,V,nset]=size(X); t_obs=zeros(1,V,nset); p_fwer=zeros(1,V,nset);
rng(seed,'twister');
for j=1:nset
    A=X(:,:,j);
    t_obs(:,:,j)=mean(A,1)./max(std(A,0,1)/sqrt(n),eps);
    mx=zeros(n_perm,1);
    for p=1:n_perm
        sg=2*(rand(n,1)>.5)-1; P=A.*sg;
        tp=mean(P,1)./max(std(P,0,1)/sqrt(n),eps);
        mx(p)=max(abs(tp));
    end
    for first=1:5000:V
        ii=first:min(first+4999,V);
        p_fwer(1,ii,j)=(1+sum(mx>=abs(t_obs(1,ii,j)),1))/(1+n_perm);
    end
end
end
