%% s018_make_summary_figures_and_tables
clear; clc;
cfg=load_project_config();
ridge_dir=cfg.result_dirs{5}; netridge_dir=cfg.result_dirs{13};
lag_dir=cfg.result_dirs{14}; out_dir=cfg.result_dirs{19};
n=numel(cfg.subject_list);
rows={};

for s=1:n
    subject_id=cfg.subject_list{s};
    M=load(fullfile(ridge_dir,sprintf( ...
        'sub-%s_ridge_main_results.mat',subject_id)), ...
        'R2_cv_voxel','prediction_correlation');
    rows(end+1,:)={subject_id,'ridge_400_zero_lag', ...
        mean(M.R2_cv_voxel),median(M.R2_cv_voxel), ...
        mean(M.R2_cv_voxel>0),mean(M.prediction_correlation)}; %#ok<SAGROW>
    for scheme=17
        f=fullfile(netridge_dir,sprintf( ...
            'sub-%s_network%d_ridge_results.mat',subject_id,scheme));
        if exist(f,'file')
            D=load(f,'result');
            rows(end+1,:)={subject_id,sprintf('ridge_network%d',scheme), ...
                mean(D.result.R2),median(D.result.R2), ...
                mean(D.result.R2>0),mean(D.result.prediction_correlation)}; %#ok<SAGROW>
        end
    end
    f=fullfile(lag_dir,sprintf('sub-%s_lagged_ridge_results.mat',subject_id));
    if exist(f,'file')
        D=load(f,'result');
        rows(end+1,:)={subject_id,'ridge_lagged',mean(D.result.R2), ...
            median(D.result.R2),mean(D.result.R2>0), ...
            mean(D.result.prediction_correlation)}; %#ok<SAGROW>
    end
end
T=cell2table(rows,'VariableNames',{'subject_id','model','mean_R2', ...
    'median_R2','positive_R2_fraction','mean_prediction_correlation'});
writetable(T,fullfile(out_dir,'model_comparison_table.csv'));

models=unique(string(T.model),'stable');
summary=table();
summary.model=models;
summary.mean_of_subject_mean_R2=zeros(numel(models),1);
summary.sd_of_subject_mean_R2=zeros(numel(models),1);
for i=1:numel(models)
    x=T.mean_R2(string(T.model)==models(i));
    summary.mean_of_subject_mean_R2(i)=mean(x);
    summary.sd_of_subject_mean_R2(i)=std(x);
end
writetable(summary,fullfile(out_dir,'model_comparison_summary.csv'));

fig=figure('Color','w','Visible','off');
boxchart(categorical(string(T.model)),T.mean_R2);
ylabel('Mean voxelwise cross-validated R^2');
title('Model comparison across subjects');
grid on;
exportgraphics(fig,fullfile(out_dir,'model_comparison.png'),'Resolution',200);
close(fig);
