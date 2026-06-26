%% FigureS01 - Mask and input QC
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS01');
fig_atlas_wm_panel(cfg.schaefer400_file,cfg.wm_mask_final, ...
    fullfile(paneldir,'panel_a_GM400_WM95_atlas.png'), ...
    'Schaefer-400 GM atlas + WM95 mask');
fig_brain_panel(cfg.wm_mask_final,fullfile(paneldir,'panel_b_group_WM95_mask.png'),'Strict WM95 mask','Reds',kind='roi');
Q=readtable(fullfile(cfg.result_dirs{2},'QC_mask_alignment_report.csv'),'Delimiter',',','TextType','string');
qn=Q; qn.value_num=str2double(string(qn.value)); qn=qn(isfinite(qn.value_num),:);
keep=ismember(qn.metric,["wm_voxels","gm_voxels","atlas_voxels","existing_runs","atlas_wm_overlap","gm_wm_overlap"]);
qn=qn(keep,:);
label=string(qn.metric);
label=strrep(label,'wm_voxels','WM voxels');
label=strrep(label,'gm_voxels','GM voxels');
label=strrep(label,'atlas_voxels','Atlas voxels');
label=strrep(label,'existing_runs','Runs');
label=strrep(label,'atlas_wm_overlap','Atlas-WM overlap');
label=strrep(label,'gm_wm_overlap','GM-WM overlap');
fig=figure('Position',[100 100 450 270]); hold on;
y=1:height(qn); x=log10(qn.value_num+1);
barh(y,x,'FaceColor',[.35 .35 .35],'EdgeColor','none');
yticks(y); yticklabels(label); xlabel('log10(value + 1)');
for i=1:height(qn)
    text(x(i)+0.05,y(i),sprintf('%g',qn.value_num(i)),'FontSize',6.5,'VerticalAlignment','middle');
end
title('Mask and run QC'); set(gca,'TickDir','out','YDir','reverse'); xlim([0 max(x)+1.0]);
fig_export(fig,fullfile(paneldir,'panel_c_mask_alignment_QC.png'));
R=readtable(fullfile(cfg.result_dirs{2},'valid_subject_run_table.csv'),'Delimiter',',','TextType','string');
fig=figure('Position',[100 100 300 240]); cats=categorical({'Runs','Subjects'}); cats=reordercats(cats,{'Runs','Subjects'}); bar(cats,[height(R) numel(unique(R.subject_id))],'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylabel('Count'); title('Included data'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_d_included_subjects_runs.png'));
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S1 layout guide\n\nArrange as 2x2: GM mask, WM mask, QC metrics, included subjects/runs.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S1 legend\n\nInput and mask quality control for the primary strict-WM analysis.\n");
fprintf('FigureS01 complete\n');
