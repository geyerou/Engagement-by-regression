function Xnet = aggregate_network_timeseries(X,network_id)
n=max(network_id);
Xnet=zeros(size(X,1),n,'single');
for k=1:n
    Xnet(:,k)=mean(X(:,network_id==k),2);
end
end
