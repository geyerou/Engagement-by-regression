function [Bz,muX,sdX,muY,sdY] = fit_final_beta_blocked( ...
    X,Y,lambda,block_size)
muX = mean(X,1);
sdX = std(X,0,1); sdX(sdX<eps('single')) = 1;
Xz = double((X-muX)./sdX);
[U,S,V] = svd(Xz,'econ');
s = diag(S);
shrink = s ./ (s.^2 + size(Xz,1)*lambda);
n_vox = size(Y,2);
Bz = zeros(size(X,2),n_vox,'single');
muY = zeros(1,n_vox,'single');
sdY = zeros(1,n_vox,'single');
for first = 1:block_size:n_vox
    idx = first:min(first+block_size-1,n_vox);
    y = Y(:,idx);
    muY(idx) = mean(y,1);
    sdY(idx) = std(y,0,1);
    sdY(idx(sdY(idx)<eps('single'))) = 1;
    yz = double((y-muY(idx))./sdY(idx));
    Bz(:,idx) = single(V*(shrink.*(U'*yz)));
end
end
