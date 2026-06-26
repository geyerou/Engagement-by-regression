function result = crossval_nested_ridge_r2(X_runs,Y_runs,lambda_grid, ...
    selection_voxels,block_size)
n_runs=numel(X_runs);
n_vox=size(Y_runs{1},2);
sample_n=min(selection_voxels,n_vox);
voxel_sample=round(linspace(1,n_vox,sample_n));
lambdas=zeros(n_runs,1);
performance=zeros(n_runs,numel(lambda_grid));
for test_run=1:n_runs
    train_runs=setdiff(1:n_runs,test_run);
    [lambdas(test_run),performance(test_run,:)]=select_lambda_inner_cv( ...
        X_runs,Y_runs,train_runs,lambda_grid,voxel_sample);
end
[R2,r]=crossval_fixed_lambda_r2(X_runs,Y_runs,lambdas,block_size);
result=struct('R2',R2,'prediction_correlation',r, ...
    'best_lambda_outer',lambdas,'lambda_performance',performance);
end
