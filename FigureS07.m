%% FigureS07 - Cross-session validation
clear; clc; code_dir=fileparts(mfilename('fullpath')); addpath(fullfile(code_dir,'functions')); fig_style(); cfg=load_project_config(); [figdir,paneldir]=fig_dir('FigureS07');
R=readtable(fullfile(cfg.output_dir,'result_s028_session_local_r2_fc_validation','session_local_R2_residual_group_summary.csv'),'Delimiter',',','TextType','string');
G=readtable(fullfile(cfg.output_dir,'result_s029_cross_session_gradient_validation','cross_session_gradient_group_summary.csv'),'Delimiter',',','TextType','string');
fig=figure('Position',[100 100 320 250]); bar(categorical({'Subject mean','Group map'}),[R.value(R.metric=="mean_subject_session_r") R.value(R.metric=="group_mean_map_r")],'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylim([0 1]); ylabel('Correlation'); title('R² residual REST1 vs REST2'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_a_session_R2_residual_validation.png'));
fig=figure('Position',[100 100 320 250]); bar(categorical({'Subject mean','Group map'}),[G.value(G.metric=="mean_subject_cross_direction_r") G.value(G.metric=="group_cross_direction_map_r")],'FaceColor',[.55 .55 .55],'EdgeColor','none'); ylim([0 1]); ylabel('Correlation'); title('Beta-FC G1 cross-session'); set(gca,'TickDir','out'); fig_export(fig,fullfile(paneldir,'panel_b_cross_session_G1_validation.png'));
fig_brain_panel(fullfile(cfg.output_dir,'result_s029_cross_session_gradient_validation','cross_session_mean_beta_minus_FC_G1.nii.gz'),fullfile(paneldir,'panel_c_cross_session_G1_map.png'),'Cross-session beta-FC G1','coolwarm',symmetric=true);
fig_write_text(fullfile(figdir,'layout_guide.md'),"# Figure S7 layout guide\n\nArrange R2 residual validation, G1 validation, and cross-session G1 map.\n");
fig_write_text(fullfile(figdir,'legend.md'),"# Figure S7 legend\n\nSession-level validation of local R2 residuals and cross-session projection validation of beta-FC G1.\n");
fprintf('FigureS07 complete\n');
