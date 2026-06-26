%% s034_g1_binned_cortical_contribution
% Bin WM voxels by the beta-FC G1 difference axis and summarize local
% prediction, FC, cortical-network contribution, and JHU composition.
clear; clc; [cfg,ext]=round1_local_extension_config();
out_dir=ext.result_dirs{14};

L=load(fullfile(ext.result_dirs{1},'local_r2_fc_group_data.mat'), ...
    'group_mean_zR2','group_mean_zFC_rms','group_mean_residual', ...
    'R2_residual','wm_idx');
G=load(fullfile(ext.result_dirs{4},'gradient_subject_inference.mat'), ...
    'mean_difference');
K=load(fullfile(ext.result_dirs{12},'spatial_confound_control.mat'), ...
    'group_R2_residual_controlled','group_G1diff_controlled', ...
    'group_R2_minus_FC_confounds');
C=load(fullfile(ext.result_dirs{11},'spatial_confound_maps.mat'), ...
    'confounds');
if isfile(fullfile(ext.result_dirs{9},'cross_session_gradient_validation.mat'))
    XG=load(fullfile(ext.result_dirs{9},'cross_session_gradient_validation.mat'), ...
        'crossfit_mean');
    bin_axis=double(XG.crossfit_mean(:));
    bin_axis_name="cross_session_beta_minus_FC_G1";
else
    bin_axis=double(G.mean_difference(:));
    bin_axis_name="subject_mean_beta_minus_FC_G1";
end

n_bins=ext.g1_bin_count;
[bin_id,edges]=equal_count_bins(bin_axis,n_bins);

% Network contribution by bin.
labels=read_schaefer_labels(cfg);
net_names=string(labels.network17_names(:));
n_net=numel(net_names);
net_sum=zeros(n_bins,n_net);
net_count=zeros(n_bins,n_net);
for s=1:numel(cfg.subject_list)
    sid=cfg.subject_list{s};
    f=fullfile(cfg.result_dirs{9},sprintf( ...
        'sub-%s_network_contribution_17net.mat',sid));
    if ~exist(f,'file')
        warning('Missing network contribution: %s',f); continue;
    end
    N=load(f,'Delta_R2','wm_voxel_indices','network_names');
    if ~isequal(N.wm_voxel_indices,L.wm_idx)
        error('WM mismatch in %s.',f);
    end
    for b=1:n_bins
        vox=bin_id==b;
        net_sum(b,:)=net_sum(b,:)+mean(double(N.Delta_R2(vox,:)),1,'omitnan');
        net_count(b,:)=net_count(b,:)+1;
    end
end
net_mean=net_sum ./ max(net_count,1);
net_fraction=net_mean ./ max(sum(abs(net_mean),2),eps);

% JHU direct-overlap composition by bin.
jhu=double(niftiread(ext.jhu_file)); jhu=jhu(L.wm_idx);
jhu_labels=read_jhu_label_names(ext.jhu_xml);
jhu_rows={};
for b=1:n_bins
    vox=bin_id==b;
    ids=unique(jhu(vox)); ids(ids==0)=[];
    for ii=1:numel(ids)
        id=ids(ii);
        name=jhu_labels.label_name(jhu_labels.label_id==id);
        if isempty(name), name="unknown"; end
        n_direct=nnz(jhu(vox)==id);
        jhu_rows(end+1,:)={b,id,name,n_direct,n_direct/nnz(vox)}; %#ok<SAGROW>
    end
end
J=cell2table(jhu_rows,'VariableNames',{'bin','label_id','label_name', ...
    'direct_voxels','fraction_within_bin'});

% Bin-level local metrics.
rows=cell(n_bins,15);
for b=1:n_bins
    vox=bin_id==b;
    [~,top_net]=max(net_fraction(b,:));
    rows(b,:)={b,nnz(vox),mean(bin_axis(vox),'omitnan'), ...
        min(bin_axis(vox)),max(bin_axis(vox)), ...
        mean(double(L.group_mean_zR2(vox)),'omitnan'), ...
        mean(double(L.group_mean_zFC_rms(vox)),'omitnan'), ...
        mean(double(L.group_mean_residual(vox)),'omitnan'), ...
        mean(double(K.group_R2_residual_controlled(vox)),'omitnan'), ...
        mean(double(K.group_R2_minus_FC_confounds(vox)),'omitnan'), ...
        mean(double(C.confounds.dist_GM_mm(vox)),'omitnan'), ...
        mean(double(C.confounds.dist_WM_boundary_mm(vox)),'omitnan'), ...
        mean(double(C.confounds.WM_prevalence(vox)),'omitnan'), ...
        net_names(top_net),net_fraction(b,top_net)};
end
B=cell2table(rows,'VariableNames',{'bin','n_voxels','mean_G1_axis', ...
    'min_G1_axis','max_G1_axis','mean_zR2','mean_zFC_RMS', ...
    'mean_R2_residual','mean_R2_residual_controlled', ...
    'mean_zR2_after_zFC_spatial_confounds','mean_dist_GM_mm', ...
    'mean_dist_WM_boundary_mm','mean_WM_prevalence', ...
    'top_network_fraction_name','top_network_fraction'});

% Long network table.
net_rows=cell(n_bins*n_net,5); rr=0;
for b=1:n_bins
    for k=1:n_net
        rr=rr+1;
        net_rows(rr,:)={b,k,net_names(k),net_mean(b,k),net_fraction(b,k)};
    end
end
NT=cell2table(net_rows,'VariableNames',{'bin','network_id', ...
    'network_name','mean_DeltaR2','fraction_abs_DeltaR2'});

writetable(B,fullfile(out_dir,'G1_bin_summary.csv'));
writetable(NT,fullfile(out_dir,'G1_bin_network_contribution.csv'));
writetable(J,fullfile(out_dir,'G1_bin_JHU_composition.csv'));

fig=figure('Color','w','Visible','off','Position',[50 50 1500 820]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; plot(B.bin,B.mean_zR2,'-o'); hold on; plot(B.bin,B.mean_zFC_RMS,'-o');
legend({'zR2','zFC RMS'},'Location','best'); xlabel('G1 bin low to high');
title('Local prediction and marginal FC');
nexttile; plot(B.bin,B.mean_R2_residual,'-o'); hold on;
plot(B.bin,B.mean_R2_residual_controlled,'-o');
plot(B.bin,B.mean_zR2_after_zFC_spatial_confounds,'-o');
legend({'R2 residual','controlled residual','zR2|zFC+conf'},'Location','best');
xlabel('G1 bin'); title('FC-controlled local prediction');
nexttile; imagesc(net_fraction'); colorbar; colormap(parula);
yticks(1:n_net); yticklabels(net_names); xlabel('G1 bin');
title('17-network fractional contribution');
nexttile; plot(B.bin,B.mean_dist_GM_mm,'-o'); hold on;
plot(B.bin,B.mean_dist_WM_boundary_mm,'-o');
legend({'GM distance','WM boundary distance'},'Location','best');
xlabel('G1 bin'); ylabel('mm'); title('Anatomical proximity');
nexttile; bar(B.bin,B.mean_WM_prevalence); ylim([0 1]);
xlabel('G1 bin'); title('Subject WM prevalence');
nexttile;
GJ=groupsummary(J,'label_name','sum','direct_voxels');
[~,ord]=sort(GJ.sum_direct_voxels,'descend');
GJ=GJ(ord,:); nshow=min(12,height(GJ));
barh(categorical(GJ.label_name(1:nshow)),GJ.sum_direct_voxels(1:nshow));
set(gca,'YDir','reverse'); xlabel('Direct WM voxels');
title('Most represented JHU labels');
exportgraphics(fig,fullfile(out_dir,'G1_bin_summary.png'),'Resolution',200);
close(fig);

save(fullfile(out_dir,'G1_binned_cortical_contribution.mat'), ...
    'B','NT','J','bin_id','edges','bin_axis','bin_axis_name','net_mean', ...
    'net_fraction','net_names','-v7.3');
fprintf('s034 complete: %d G1 bins saved to %s\n',n_bins,out_dir);

function [bin_id,edges]=equal_count_bins(x,n_bins)
x=x(:); [xs,ord]=sort(x);
bin_id=zeros(size(x));
edges=zeros(n_bins+1,1); edges(1)=xs(1); edges(end)=xs(end);
for b=1:n_bins
    lo=floor((b-1)*numel(x)/n_bins)+1;
    hi=floor(b*numel(x)/n_bins);
    idx=ord(lo:hi);
    bin_id(idx)=b;
    edges(b)=xs(lo); edges(b+1)=xs(hi);
end
end
