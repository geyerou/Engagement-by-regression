function [R2,prediction_correlation] = crossval_fixed_lambda_r2( ...
    X_runs,Y_runs,lambdas,block_size)
n_runs = numel(X_runs);
n_vox = size(Y_runs{1},2);
sse=zeros(1,n_vox); sst=zeros(1,n_vox);
sy=zeros(1,n_vox); sp=zeros(1,n_vox);
syy=zeros(1,n_vox); spp=zeros(1,n_vox); syp=zeros(1,n_vox);
n_total=0;
for test_run=1:n_runs
    train_runs=setdiff(1:n_runs,test_run);
    pred=ridge_fit_predict_blocked(vertcat(X_runs{train_runs}), ...
        vertcat(Y_runs{train_runs}),X_runs{test_run}, ...
        lambdas(test_run),block_size,false);
    y=double(Y_runs{test_run}); p=double(pred);
    sse=sse+sum((y-p).^2,1);
    sst=sst+sum((y-mean(y,1)).^2,1);
    sy=sy+sum(y,1); sp=sp+sum(p,1);
    syy=syy+sum(y.^2,1); spp=spp+sum(p.^2,1); syp=syp+sum(y.*p,1);
    n_total=n_total+size(y,1);
end
R2=single(1-sse./max(sst,eps));
num=n_total*syp-sy.*sp;
den=sqrt(max(n_total*syy-sy.^2,0).*max(n_total*spp-sp.^2,0));
prediction_correlation=single(num./max(den,eps));
end
