function [prediction,beta] = ridge_fit_predict_blocked( ...
    Xtrain,Ytrain,Xtest,lambda,block_size,return_beta)
if nargin < 7, return_beta = false; end
n_vox = size(Ytrain,2);
prediction = zeros(size(Xtest,1),n_vox,'single');
if return_beta
    beta = zeros(size(Xtrain,2),n_vox,'single');
else
    beta = [];
end

muX = mean(Xtrain,1);
sdX = std(Xtrain,0,1); sdX(sdX<eps('single')) = 1;
Xtr = double((Xtrain-muX)./sdX);
Xte = double((Xtest-muX)./sdX);
[U,S,V] = svd(Xtr,'econ');
s = diag(S);
shrink = s ./ (s.^2 + size(Xtr,1)*lambda);

for first = 1:block_size:n_vox
    idx = first:min(first+block_size-1,n_vox);
    Y = Ytrain(:,idx);
    muY = mean(Y,1);
    sdY = std(Y,0,1); sdY(sdY<eps('single')) = 1;
    Yz = double((Y-muY)./sdY);
    Bz = V * (shrink .* (U' * Yz));
    pred = (Xte*Bz).*double(sdY) + double(muY);
    prediction(:,idx) = single(pred);
    if return_beta
        % Standardized beta is the comparable predictive fingerprint.
        beta(:,idx) = single(Bz);
    end
end
end
