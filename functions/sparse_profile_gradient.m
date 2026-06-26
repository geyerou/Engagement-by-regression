function [gradient,eigenvalues] = sparse_profile_gradient( ...
    profile,n_components,n_pca,k_neighbors)
% Scalable cosine/normalized-angle kNN diffusion embedding.
X=double(profile);
X=X-mean(X,2);
X=X./max(sqrt(sum(X.^2,2)),eps);
n_pca=min([n_pca,size(X,2),size(X,1)-1]);
[~,score]=pca(X,'NumComponents',n_pca,'Centered',true);
score=score./max(sqrt(sum(score.^2,2)),eps);
[idx,dist]=knnsearch(score,score,'K',k_neighbors+1, ...
    'Distance','cosine');
idx=idx(:,2:end); dist=dist(:,2:end);
n=size(X,1); row=repelem((1:n)',k_neighbors,1); col=idx(:);
cosine=max(-1,min(1,1-dist(:)));
weight=1-acos(cosine)/pi;
W=sparse(row,col,weight,n,n); W=max(W,W');
d=sum(W,2); A=spdiags(1./sqrt(max(d,eps)),0,n,n)*W* ...
    spdiags(1./sqrt(max(d,eps)),0,n,n);
[V,D]=eigs(A,n_components+1,'largestreal', ...
    struct('tol',1e-5,'maxit',1000));
[ev,ord]=sort(diag(D),'descend'); V=V(:,ord);
gradient=single(V(:,2:n_components+1));
eigenvalues=double(ev(2:n_components+1));
end
