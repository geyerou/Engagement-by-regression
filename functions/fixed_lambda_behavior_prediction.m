function result = fixed_lambda_behavior_prediction( ...
    X,C,y,outer_folds,selected_lambda)
% Repeated out-of-sample behavior prediction using preselected penalties.

n = numel(y);
n_repeats = size(outer_folds,2);
pred_full_repeats = zeros(n,n_repeats);
pred_confound_repeats = zeros(n,n_repeats);
lambda_counter = 0;

for r = 1:n_repeats
    fold_labels = outer_folds(:,r);
    for outer = unique(fold_labels)'
        lambda_counter = lambda_counter+1;
        lambda = selected_lambda(lambda_counter);
        test = fold_labels==outer;
        train = ~test;

        Xtr = X(train,:); Xte = X(test,:);
        Ctr = C(train,:); Cte = C(test,:);
        ytr = y(train);

        b_y = Ctr\ytr;
        pred_confound_repeats(test,r) = Cte*b_y;
        yres = ytr-Ctr*b_y;

        b_x = Ctr\Xtr;
        Ztr = Xtr-Ctr*b_x;
        Zte = Xte-Cte*b_x;
        mu = mean(Ztr,1);
        sd = std(Ztr,0,1); sd(sd<eps) = 1;
        Ztr = (Ztr-mu)./sd;
        Zte = (Zte-mu)./sd;

        [U,S,V] = svd(Ztr,'econ');
        sv = diag(S);
        shrink = sv./(sv.^2+size(Ztr,1)*lambda);
        pred_full_repeats(test,r) = Cte*b_y + ...
            Zte*(V*(shrink.*(U'*yres)));
    end
end
if lambda_counter ~= numel(selected_lambda)
    error('Selected-lambda count does not match outer folds.');
end

pred_full = mean(pred_full_repeats,2);
pred_confound = mean(pred_confound_repeats,2);
sst = sum((y-mean(y)).^2);
result = struct();
result.full_r = corr(y,pred_full);
result.confound_r = corr(y,pred_confound);
result.full_R2 = 1-sum((y-pred_full).^2)/max(sst,eps);
result.confound_R2 = 1-sum((y-pred_confound).^2)/max(sst,eps);
result.delta_r = result.full_r-result.confound_r;
result.delta_R2 = result.full_R2-result.confound_R2;
end
