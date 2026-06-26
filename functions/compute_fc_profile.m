function FC = compute_fc_profile(X,Y,block_size)
% Compute voxel x predictor Pearson correlations in memory-safe blocks.
X = single(X);
Y = single(Y);
X = X-mean(X,1);
X = X./max(std(X,0,1),eps('single'));
n_vox = size(Y,2);
FC = zeros(n_vox,size(X,2),'single');
for first = 1:block_size:n_vox
    idx = first:min(first+block_size-1,n_vox);
    y = Y(:,idx);
    y = y-mean(y,1);
    y = y./max(std(y,0,1),eps('single'));
    FC(idx,:) = single((double(y)'*double(X))/(size(X,1)-1));
end
end
