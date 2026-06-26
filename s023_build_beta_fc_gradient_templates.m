%% s023_build_beta_fc_gradient_templates
% Fixed Round-1 beta/FC G1 templates and linear fingerprint projections.
% The templates are independent split averages from completed s015b.
clear; clc; [cfg,ext]=round1_local_extension_config();
wm=niftiread(cfg.wm_mask_final)>0; wm_idx=find(wm);
gdir=fullfile(cfg.result_dirs{16},'s015b_gradient_stability');

ba=niftiread(fullfile(gdir,'session_beta_splitA_aligned_gradient-01.nii.gz'));
bb=niftiread(fullfile(gdir,'session_beta_splitB_aligned_gradient-01.nii.gz'));
fa=niftiread(fullfile(gdir,'session_FC_splitA_aligned_gradient-01.nii.gz'));
fb=niftiread(fullfile(gdir,'session_FC_splitB_aligned_gradient-01.nii.gz'));
beta_template=mean([double(ba(wm_idx)),double(bb(wm_idx))],2);
fc_template=mean([double(fa(wm_idx)),double(fb(wm_idx))],2);
if corr(beta_template,fc_template)<0,fc_template=-fc_template;end
beta_template=zscore(beta_template); fc_template=zscore(fc_template);
template_difference=single(beta_template-fc_template);

G=load(fullfile(cfg.result_dirs{12},'s011b_fc_vs_beta', ...
    'FC_vs_beta_group_results.mat'),'group_mean_beta','group_mean_fc');
[beta_projection,beta_mu,beta_sd,beta_fit_r]=fit_projection( ...
    double(G.group_mean_beta),beta_template,ext.gradient_projection_lambda);
[fc_projection,fc_mu,fc_sd,fc_fit_r]=fit_projection( ...
    double(G.group_mean_fc),fc_template,ext.gradient_projection_lambda);

write_map_like(beta_template,wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{3},'beta_G1_fixed_template.nii.gz'));
write_map_like(fc_template,wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{3},'FC_G1_fixed_template.nii.gz'));
write_map_like(template_difference,wm_idx,cfg.wm_mask_final, ...
    fullfile(ext.result_dirs{3},'beta_minus_FC_G1_template.nii.gz'));
metric=["beta_split_r";"FC_split_r";"beta_FC_template_r"; ...
    "beta_projection_fit_r";"FC_projection_fit_r"];
value=[corr(double(ba(wm_idx)),double(bb(wm_idx))); ...
    corr(double(fa(wm_idx)),double(fb(wm_idx))); ...
    corr(beta_template,fc_template);beta_fit_r;fc_fit_r];
writetable(table(metric,value),fullfile(ext.result_dirs{3}, ...
    'gradient_template_summary.csv'));
save(fullfile(ext.result_dirs{3},'gradient_templates_and_projections.mat'), ...
    'beta_template','fc_template','template_difference', ...
    'beta_projection','fc_projection','beta_mu','beta_sd','fc_mu','fc_sd', ...
    'beta_fit_r','fc_fit_r','wm_idx','-v7.3');
fprintf('s023 complete: template r=%.3f; projection fits beta %.3f FC %.3f\n', ...
    corr(beta_template,fc_template),beta_fit_r,fc_fit_r);

function [b,mu,sd,fit_r]=fit_projection(X,y,lambda)
mu=mean(X,1); sd=std(X,0,1); sd(sd<eps)=1;
Z=(X-mu)./sd; b=(Z'*Z+size(Z,1)*lambda*eye(size(Z,2)))\(Z'*y);
fit_r=corr(Z*b,y);
end
