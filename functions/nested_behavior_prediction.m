function result = nested_behavior_prediction(X,C,y,outer_folds, ...
    inner_k,lambda_grid,seed)
% Repeated nested CV with training-fold confound adjustment.
%
% Confounds are unpenalized. Imaging features are residualized against
% confounds and standardized using training data only.

n = numel(y);
n_repeats = size(outer_folds,2);
pred_full_repeats = zeros(n,n_repeats);
pred_confound_repeats = zeros(n,n_repeats);
selected_lambda = [];

for r = 1:n_repeats
    fold_labels = outer_folds(:,r);
    unique_folds = unique(fold_labels)';
    for outer = unique_folds
        test = fold_labels==outer;
        train = ~test;

        Xtr = X(train,:); Xte = X(test,:);
        Ctr = C(train,:); Cte = C(test,:);
        ytr = y(train); yte = y(test); %#ok<NASGU>

        % Confound-only prediction.
        b_c = Ctr\ytr;
        pred_confound_repeats(test,r) = Cte*b_c;

        % Inner folds are stratified approximately by quantiles of y and
        % generated only within the outer training set.
        y_bins = discretize(ytr, ...
            unique(quantile(ytr,linspace(0,1,min(inner_k,5)+1))));
        if any(isnan(y_bins)) || numel(unique(y_bins))<2
            y_bins = ones(size(ytr));
        end
        inner_labels = make_repeated_stratified_folds( ...
            string(y_bins),inner_k,1,seed+100*r+outer);
        inner_labels = inner_labels(:,1);
        mse = zeros(numel(lambda_grid),1);
        count = zeros(numel(lambda_grid),1);

        for inner = unique(inner_labels)'
            iva = inner_labels==inner;
            itr = ~iva;
            [Ztr,Zva,yres_tr,conf_va] = prepare_behavior_fold( ...
                Xtr(itr,:),Xtr(iva,:),Ctr(itr,:),Ctr(iva,:), ...
                ytr(itr,:));
            [U,S,V] = svd(Ztr,'econ');
            sv = diag(S);
            Uy = U'*yres_tr;
            for l = 1:numel(lambda_grid)
                shrink = sv./(sv.^2+size(Ztr,1)*lambda_grid(l));
                pred = conf_va + Zva*(V*(shrink.*Uy));
                err = ytr(iva)-pred;
                mse(l) = mse(l)+sum(err.^2);
                count(l) = count(l)+numel(err);
            end
        end
        mse = mse./max(count,1);
        [~,best] = min(mse);
        lambda = lambda_grid(best);
        selected_lambda(end+1,1) = lambda; %#ok<AGROW>

        [Ztr,Zte,yres_tr,conf_te] = prepare_behavior_fold( ...
            Xtr,Xte,Ctr,Cte,ytr);
        [U,S,V] = svd(Ztr,'econ');
        sv = diag(S);
        shrink = sv./(sv.^2+size(Ztr,1)*lambda);
        pred_full_repeats(test,r) = conf_te + ...
            Zte*(V*(shrink.*(U'*yres_tr)));
    end
end

pred_full = mean(pred_full_repeats,2);
pred_confound = mean(pred_confound_repeats,2);
sst = sum((y-mean(y)).^2);

result = struct();
result.pred_full = pred_full;
result.pred_confound = pred_confound;
result.pred_full_repeats = pred_full_repeats;
result.pred_confound_repeats = pred_confound_repeats;
result.full_r = corr(y,pred_full);
result.confound_r = corr(y,pred_confound);
result.full_R2 = 1-sum((y-pred_full).^2)/max(sst,eps);
result.confound_R2 = 1-sum((y-pred_confound).^2)/max(sst,eps);
result.full_MAE = mean(abs(y-pred_full));
result.confound_MAE = mean(abs(y-pred_confound));
result.delta_r = result.full_r-result.confound_r;
result.delta_R2 = result.full_R2-result.confound_R2;
result.selected_lambda = selected_lambda;
end

function [Ztr,Zte,yres_tr,conf_te] = prepare_behavior_fold( ...
    Xtr,Xte,Ctr,Cte,ytr)
% Training-only confound regression and standardization.
b_y = Ctr\ytr;
yres_tr = ytr-Ctr*b_y;
conf_te = Cte*b_y;

b_x = Ctr\Xtr;
Ztr = Xtr-Ctr*b_x;
Zte = Xte-Cte*b_x;
mu = mean(Ztr,1);
sd = std(Ztr,0,1);
sd(sd<eps) = 1;
Ztr = (Ztr-mu)./sd;
Zte = (Zte-mu)./sd;
end
