function [best_lambda,mean_score] = select_lambda_inner_cv( ...
    X_runs,Y_runs,train_run_indices,lambda_grid,voxel_sample)
% Select one subject/fold-level lambda using only outer-training runs.
n_lambda = numel(lambda_grid);
scores = zeros(numel(train_run_indices),n_lambda);

for ii = 1:numel(train_run_indices)
    val_run = train_run_indices(ii);
    inner_train = setdiff(train_run_indices,val_run);
    Xtr = vertcat(X_runs{inner_train});
    Xva = X_runs{val_run};
    Ytr = vertcat(Y_runs{inner_train});
    Yva = Y_runs{val_run};
    Ytr = Ytr(:,voxel_sample);
    Yva = Yva(:,voxel_sample);
    [XtrZ,XvaZ,YtrZ,~,muY,sdY] = ...
        standardize_train_test(Xtr,Xva,Ytr,Yva);
    [U,S,V] = svd(double(XtrZ),'econ');
    singular_values = diag(S);
    UY = U' * double(YtrZ);
    sst = sum((double(Yva)-mean(double(Yva),1)).^2,'all');
    for l = 1:n_lambda
        shrink = singular_values ./ ...
            (singular_values.^2 + size(XtrZ,1)*lambda_grid(l));
        predZ = double(XvaZ) * (V * (shrink .* UY));
        pred = predZ .* double(sdY) + double(muY);
        sse = sum((double(Yva)-pred).^2,'all');
        scores(ii,l) = 1-sse/max(sst,eps);
    end
end
mean_score = mean(scores,1);
[~,idx] = max(mean_score);
best_lambda = lambda_grid(idx);
end
