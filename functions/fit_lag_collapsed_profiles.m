function [integrated_beta,lag_energy,lag_center,best_lag] = ...
    fit_lag_collapsed_profiles(X,Y,lambda,block_size,n_roi,lags)
% Fit full-data standardized lag model but retain only compact voxel profiles.
muX=mean(X,1); sdX=std(X,0,1); sdX(sdX<eps('single'))=1;
Xz=double((X-muX)./sdX); [U,S,V]=svd(Xz,'econ');
sv=diag(S); shrink=sv./(sv.^2+size(Xz,1)*lambda);
n_vox=size(Y,2); L=numel(lags);
integrated_beta=zeros(n_vox,n_roi,'single');
lag_energy=zeros(n_vox,L,'single');
lag_center=zeros(n_vox,1,'single');
best_lag=zeros(n_vox,1,'single');
for first=1:block_size:n_vox
    idx=first:min(first+block_size-1,n_vox);
    y=Y(:,idx); my=mean(y,1); sy=std(y,0,1);
    sy(sy<eps('single'))=1;
    B=single(V*(shrink.*(U'*double((y-my)./sy))));
    B3=reshape(B,n_roi,L,numel(idx));
    integrated_beta(idx,:)=squeeze(sum(B3,2))';
    E=squeeze(sqrt(sum(B3.^2,1)))';
    lag_energy(idx,:)=E;
    w=E.^2; lag_center(idx)=single((w*double(lags(:)))./max(sum(w,2),eps));
    [~,ii]=max(E,[],2); best_lag(idx)=single(lags(ii));
end
end
