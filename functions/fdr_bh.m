function q = fdr_bh(p)
shape=size(p); p=p(:);
[ps,order]=sort(p);
m=numel(p);
qs=ps.*m./(1:m)';
qs=flipud(cummin(flipud(qs)));
qs=min(qs,1);
q=zeros(m,1); q(order)=qs;
q=reshape(q,shape);
end
