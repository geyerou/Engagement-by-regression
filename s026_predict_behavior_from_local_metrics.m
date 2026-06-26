%% s026_predict_behavior_from_local_metrics
% Nested prediction with fixed anatomical features and Freedman-Lane tests.
clear; clc; [cfg,ext]=round1_local_extension_config();
rng(20260620,'twister');
wm=niftiread(cfg.wm_mask_final)>0; wm_idx=find(wm);
jhu=double(niftiread(ext.jhu_file)); labels=read_jhu_label_names(ext.jhu_xml);
counts=arrayfun(@(k)nnz(jhu(wm_idx)==k),labels.label_id);
tract_ids=labels.label_id(counts>=20);
tract_names=labels.label_name(counts>=20);
n=numel(cfg.subject_list); subject_id=string(cfg.subject_list(:));
Xr=zeros(n,numel(tract_ids)); Xg=zeros(n,numel(tract_ids));
for s=1:n
    sid=cfg.subject_list{s};
    R=load(fullfile(ext.result_dirs{1},sprintf( ...
        'sub-%s_local_r2_fc.mat',sid)),'residual');
    G=load(fullfile(ext.result_dirs{4},sprintf( ...
        'sub-%s_beta_FC_G1_scores.mat',sid)),'zb','zf');
    gd=G.zb-G.zf;
    for k=1:numel(tract_ids)
        m=jhu(wm_idx)==tract_ids(k);
        Xr(s,k)=mean(R.residual(m),'omitnan');
        Xg(s,k)=mean(gd(m),'omitnan');
    end
end
feature_sets={Xr,Xg,[Xr Xg]};
feature_set_names=["R2_residual","beta_minus_FC_G1","combined"];

B=readtable(cfg.behavior_table,'TextType','string'); B.Subject=string(B.Subject);
M=readtable(cfg.motion_table,'TextType','string'); M.Subject=string(M.Subject);
[fb,ib]=ismember(subject_id,B.Subject); [fm,im]=ismember(subject_id,M.Subject);
if ~all(fb) || ~all(fm),error('Behavior/motion subjects missing.');end
B=B(ib,:); M=M(im,:);
age=categorical(B.Age); sex=categorical(B.Gender);
age_levels=categories(age); sex_levels=categories(sex);
C=ones(n,1);
for k=2:numel(age_levels),C(:,end+1)=double(age==age_levels{k});end %#ok<SAGROW>
for k=2:numel(sex_levels),C(:,end+1)=double(sex==sex_levels{k});end %#ok<SAGROW>
C(:,end+1)=mean(double(M{:,{'rfMRI_REST1_LR','rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR','rfMRI_REST2_RL'}}),2);
strata=string(B.Age)+"_"+string(B.Gender);
outer_folds=make_repeated_stratified_folds(strata, ...
    ext.behavior_outer_folds,ext.behavior_outer_repeats,20260620);

n_test=numel(feature_sets)*numel(ext.behavior_outcomes);
rows=cell(n_test,15); store=struct(); counter=0;
for fs=1:numel(feature_sets)
    X=feature_sets{fs};
    for o=1:numel(ext.behavior_outcomes)
        outcome=ext.behavior_outcomes{o}; y=double(B.(outcome));
        valid=isfinite(y)&all(isfinite(X),2)&all(isfinite(C),2);
        observed=nested_behavior_prediction(X(valid,:),C(valid,:),y(valid), ...
            outer_folds(valid,:),ext.behavior_inner_folds, ...
            ext.behavior_lambda_grid,20260620+100*fs+o);
        Cv=C(valid,:); yv=y(valid);
        fit=Cv*(Cv\yv); resid=yv-fit;
        null_delta_r=zeros(ext.behavior_n_permutations,1);
        null_delta_R2=zeros(ext.behavior_n_permutations,1);
        for p=1:ext.behavior_n_permutations
            yp=fit+resid(randperm(numel(resid)));
            Q=fixed_lambda_behavior_prediction(X(valid,:),Cv,yp, ...
                outer_folds(valid,:),observed.selected_lambda);
            null_delta_r(p)=Q.delta_r; null_delta_R2(p)=Q.delta_R2;
        end
        pr=(1+sum(null_delta_r>=observed.delta_r))/ ...
            (1+ext.behavior_n_permutations);
        pR2=(1+sum(null_delta_R2>=observed.delta_R2))/ ...
            (1+ext.behavior_n_permutations);
        counter=counter+1;
        rows(counter,:)={feature_set_names(fs),outcome,nnz(valid), ...
            size(X,2),observed.full_r,observed.full_R2, ...
            observed.confound_r,observed.confound_R2, ...
            observed.delta_r,observed.delta_R2,pr,pR2, ...
            median(observed.selected_lambda), ...
            prctile(observed.selected_lambda,25), ...
            prctile(observed.selected_lambda,75)};
        key=matlab.lang.makeValidName(feature_set_names(fs)+"_"+outcome);
        store.(key)=struct('observed',observed, ...
            'null_delta_r',null_delta_r,'null_delta_R2',null_delta_R2);
    end
end
T=cell2table(rows,'VariableNames',{'feature_set','outcome','N','n_features', ...
    'full_r','full_R2','confound_r','confound_R2','delta_r','delta_R2', ...
    'permutation_p_delta_r','permutation_p_delta_R2', ...
    'median_lambda','lambda_p25','lambda_p75'});
T.permutation_q_delta_r=fdr_bh(T.permutation_p_delta_r);
T.permutation_q_delta_R2=fdr_bh(T.permutation_p_delta_R2);
writetable(T,fullfile(ext.result_dirs{6}, ...
    'behavior_local_prediction_summary.csv'));
feature_info=table(tract_ids,tract_names,counts(counts>=20), ...
    'VariableNames',{'JHU_label','JHU_name','WM_voxels'});
writetable(feature_info,fullfile(ext.result_dirs{6}, ...
    'behavior_local_feature_tracts.csv'));
save(fullfile(ext.result_dirs{6},'behavior_local_prediction_results.mat'), ...
    'T','store','Xr','Xg','feature_info','subject_id','C','outer_folds','-v7.3');
fprintf('s026 complete: %d nested behavior tests.\n',height(T));
