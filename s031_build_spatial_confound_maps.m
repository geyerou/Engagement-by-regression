%% s031_build_spatial_confound_maps
% Build voxelwise nuisance/proximity maps for strict WM95 analyses.
% These maps are used to ask whether the local R2 residual and beta-FC G1
% effects survive simple anatomical/signal confounds.
clear; clc; [cfg,ext]=round1_local_extension_config();
out_dir=ext.result_dirs{11};

wm=logical(niftiread(cfg.wm_mask_final)); wm_idx=find(wm);
gm=logical(niftiread(cfg.gm_mask_source));
info=niftiinfo(cfg.wm_mask_final);
pix=double(info.PixelDimensions(1:3));
V=numel(wm_idx);

% Distance to the nearest group GM voxel, in mm.  This is the most direct
% gray-matter contamination/proximity covariate.
dist_gm_vox=bwdist(gm);
dist_gm_mm=single(dist_gm_vox(wm_idx) .* mean(pix));

% Distance to the WM boundary, in mm.  Boundary-proximal voxels are more
% vulnerable to partial volume and mask-normalized smoothing edge effects.
dist_nonwm_vox=bwdist(~wm);
dist_nonwm_mm=single(dist_nonwm_vox(wm_idx) .* mean(pix));

% Subject-level WM prevalence at each group-WM voxel.
wm_count=zeros(size(wm),'single');
for s=1:numel(cfg.subject_list)
    sid=cfg.subject_list{s};
    f=fullfile(cfg.data_dir,'Masks',sprintf('%s_WM_mask.nii.gz',sid));
    if exist(f,'file')
        wm_count=wm_count+single(niftiread(f)>0);
    else
        warning('Missing subject WM mask: %s',f);
    end
end
wm_prevalence=single(wm_count(wm_idx) ./ numel(cfg.subject_list));

% Temporal SD of the extracted/smoothed WM signal.  The extracted data were
% mean-centered, so this is not tSNR; it is a practical signal/noise proxy.
temporal_sd=zeros(V,1,'single');
run_count=0;
for s=1:numel(cfg.subject_list)
    sid=cfg.subject_list{s};
    for r=1:numel(cfg.run_ids)
        f=fullfile(cfg.result_dirs{3},sprintf( ...
            'sub-%s_run-%s_timeseries.mat',sid,cfg.run_ids{r}));
        if ~exist(f,'file'), warning('Missing extracted run: %s',f); continue; end
        R=load(f,'Y_wm_raw','wm_voxel_indices');
        if ~isequal(R.wm_voxel_indices,wm_idx)
            error('WM index mismatch in %s.',f);
        end
        temporal_sd=temporal_sd+single(std(double(R.Y_wm_raw),0,1)');
        run_count=run_count+1;
    end
end
temporal_sd=temporal_sd ./ max(run_count,1);

% Optional CSF/ventricle distance: only computed if a matching group CSF mask
% exists in the dataset.  Most current HCP mask folders do not contain one.
csf_file=find_existing_file({ ...
    fullfile(cfg.data_dir,'Group_CSF_Mask.nii.gz'), ...
    fullfile(cfg.data_dir,'Group_CSF_mask.nii.gz'), ...
    fullfile(cfg.data_dir,'Group_Ventricle_Mask.nii.gz'), ...
    fullfile(cfg.data_dir,'Masks','Group_CSF_Mask.nii.gz')});
if strlength(csf_file)>0
    csf=logical(niftiread(csf_file));
    dist_csf_vox=bwdist(csf);
    dist_csf_mm=single(dist_csf_vox(wm_idx) .* mean(pix));
else
    dist_csf_mm=single(nan(V,1));
end

% World coordinates.  If the NIfTI transform is unavailable/inconvenient,
% fall back to voxel coordinates in mm; the relative spatial trends are still
% useful for confound control.
[ii,jj,kk]=ind2sub(size(wm),wm_idx);
[x_mm,y_mm,z_mm]=voxel_world_coordinates(info,ii,jj,kk,pix);

write_map_like(dist_gm_mm,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'distance_to_group_GM_mm.nii.gz'));
write_map_like(dist_nonwm_mm,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'distance_to_WM_boundary_mm.nii.gz'));
write_map_like(wm_prevalence,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'subject_WM_mask_prevalence.nii.gz'));
write_map_like(temporal_sd,wm_idx,cfg.wm_mask_final, ...
    fullfile(out_dir,'mean_extracted_WM_temporal_SD.nii.gz'));
if any(isfinite(dist_csf_mm))
    write_map_like(dist_csf_mm,wm_idx,cfg.wm_mask_final, ...
        fullfile(out_dir,'distance_to_CSF_or_ventricle_mm.nii.gz'));
end

confounds=table(wm_idx(:),dist_gm_mm,dist_nonwm_mm,wm_prevalence, ...
    temporal_sd,dist_csf_mm,single(x_mm(:)),single(y_mm(:)),single(z_mm(:)), ...
    'VariableNames',{'wm_index','dist_GM_mm','dist_WM_boundary_mm', ...
    'WM_prevalence','temporal_SD','dist_CSF_mm','x_mm','y_mm','z_mm'});
writetable(confounds,fullfile(out_dir,'spatial_confound_voxel_table.csv'));

metric=["n_WM_voxels";"n_subjects";"n_runs_used_for_temporal_SD"; ...
    "median_dist_GM_mm";"median_dist_WM_boundary_mm"; ...
    "median_WM_prevalence";"median_temporal_SD"; ...
    "CSF_distance_available"];
value=[V;numel(cfg.subject_list);run_count;median(dist_gm_mm,'omitnan'); ...
    median(dist_nonwm_mm,'omitnan');median(wm_prevalence,'omitnan'); ...
    median(temporal_sd,'omitnan');double(any(isfinite(dist_csf_mm)))];
writetable(table(metric,value),fullfile(out_dir,'spatial_confound_summary.csv'));

save(fullfile(out_dir,'spatial_confound_maps.mat'),'confounds','wm_idx', ...
    'dist_gm_mm','dist_nonwm_mm','wm_prevalence','temporal_sd', ...
    'dist_csf_mm','x_mm','y_mm','z_mm','csf_file','-v7.3');
fprintf('s031 complete: spatial confound maps saved to %s\n',out_dir);

function f=find_existing_file(candidates)
f="";
for i=1:numel(candidates)
    if exist(candidates{i},'file'), f=string(candidates{i}); return; end
end
end

function [x,y,z]=voxel_world_coordinates(info,i,j,k,pix)
ijk=[double(i(:)) double(j(:)) double(k(:)) ones(numel(i),1)];
try
    [x,y,z]=transformPointsForward(info.Transform,double(i(:)), ...
        double(j(:)),double(k(:)));
    return;
catch
end
try
    T=info.Transform.T;
    xyz=ijk*T;
    x=xyz(:,1); y=xyz(:,2); z=xyz(:,3);
catch
    x=(double(i(:))-mean(double(i(:))))*pix(1);
    y=(double(j(:))-mean(double(j(:))))*pix(2);
    z=(double(k(:))-mean(double(k(:))))*pix(3);
end
end
