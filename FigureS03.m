%% FigureS03 - Full s014 circular-shift null
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS03');
wm=niftiread(cfg.wm_mask_final)>0; qfiles=dir(fullfile(cfg.result_dirs{15},'*_null_qmap_R2.nii.gz')); qfrac=zeros(numel(qfiles),1); pfrac=zeros(numel(qfiles),1); qcons=zeros(nnz(wm),1,'single'); idx=find(wm);
for i=1:numel(qfiles)
    q=niftiread(fullfile(qfiles(i).folder,qfiles(i).name)); sid=extractBetween(qfiles(i).name,'sub-','_null_qmap_R2');
    p=niftiread(fullfile(qfiles(i).folder,sprintf('sub-%s_null_pmap_R2.nii.gz',sid{1})));
    qfrac(i)=mean(q(wm)<.05); pfrac(i)=mean(p(wm)<.05); qcons=qcons+single(q(idx)<.05);
end
qcons=qcons/numel(qfiles); write_map_like(qcons,idx,cfg.wm_mask_final,fullfile(paneldir,'panel_c_fraction_subjects_q05_map.nii.gz'));
fig=figure('Position',[100 100 330 260]); histogram(pfrac,20,'FaceColor',[.33 .44 .68],'EdgeColor','none'); xlabel('Fraction p<0.05'); ylabel('Subjects'); title('Voxelwise null p extent'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_a_p05_extent_distribution.png'));
fig=figure('Position',[100 100 330 260]); histogram(qfrac,20,'FaceColor',[.20 .55 .50],'EdgeColor','none'); xlabel('Fraction q_null < 0.05'); ylabel('Subjects'); title('Voxelwise null FDR extent'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_b_q05_extent_distribution.png'));
fig_brain_panel(fullfile(paneldir,'panel_c_fraction_subjects_q05_map.nii.gz'),fullfile(paneldir,'panel_c_fraction_subjects_q05_map.png'),'Fraction subjects q_null < 0.05','viridis',vmax=1);
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S3 layout guide\n\nArrange p<0.05 extent, q<0.05 extent, and cross-subject q<0.05 consistency map.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S3 legend\n\nFull s014 circular-shift null validation across all 81 subjects.\n");
fprintf('FigureS03 complete\n');
