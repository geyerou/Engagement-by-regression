function [Xbasis,Ytrim,basis] = build_lagged_basis_run(X,Y,lags,n_basis)
% Smooth DCT representation of the lag axis. Constant plus low-frequency
% cosine components replace one independent coefficient per lag.
[Xlag,Ytrim]=build_lagged_run(X,Y,lags);
L=numel(lags); R=size(X,2); T=size(Xlag,1);
basis=zeros(L,n_basis);
for k=0:n_basis-1
    basis(:,k+1)=cos(pi*((0:L-1)'+0.5)*k/L);
end
basis=basis./sqrt(sum(basis.^2,1));
X3=reshape(Xlag,T,R,L);
Xbasis=zeros(T,R*n_basis,'single');
for k=1:n_basis
    Xbasis(:,(k-1)*R+(1:R))=sum(X3.* ...
        reshape(single(basis(:,k)),1,1,L),3);
end
end
