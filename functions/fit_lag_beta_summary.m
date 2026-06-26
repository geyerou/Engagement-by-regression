function [lag_energy,best_lag] = fit_lag_beta_summary( ...
    X,Y,lambda,block_size,n_features,n_lags,lags)
muX=mean(X,1); sdX=std(X,0,1); sdX(sdX<eps('single'))=1;
Xz=double((X-muX)./sdX);
[U,S,V]=svd(Xz,'econ');
s=diag(S); shrink=s./(s.^2+size(Xz,1)*lambda);
n_vox=size(Y,2);
lag_energy=zeros(n_vox,n_lags,'single');
for first=1:block_size:n_vox
    idx=first:min(first+block_size-1,n_vox);
    y=Y(:,idx); muy=mean(y,1); sdy=std(y,0,1);
    sdy(sdy<eps('single'))=1;
    yz=double((y-muy)./sdy);
    B=single(V*(shrink.*(U'*yz)));
    for j=1:n_lags
        rows=(j-1)*n_features+(1:n_features);
        lag_energy(idx,j)=sqrt(sum(B(rows,:).^2,1))';
    end
end
[~,which]=max(lag_energy,[],2);
best_lag=single(lags(which));
end
