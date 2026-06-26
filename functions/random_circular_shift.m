function shifted = random_circular_shift(Y,min_shift)
T=size(Y,1);
allowed=min_shift:(T-min_shift);
if isempty(allowed)
    error('min_shift is too large for run length.');
end
offset=allowed(randi(numel(allowed)));
shifted=circshift(Y,offset,1);
end
