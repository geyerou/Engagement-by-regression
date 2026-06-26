function integrated_beta = fit_basis_integrated_profile( ...
    X,Y,lambda,block_size,n_roi,basis)
% Collapse smooth-basis coefficients into the integral over reconstructed lags.
[B,~,~,~,~]=fit_final_beta_blocked(X,Y,lambda,block_size);
K=size(basis,2); n_vox=size(Y,2);
B3=reshape(B,n_roi,K,n_vox);
w=single(sum(basis,1));
integrated_beta=zeros(n_vox,n_roi,'single');
for k=1:K
    integrated_beta=integrated_beta+squeeze(B3(:,k,:))'.*w(k);
end
end
