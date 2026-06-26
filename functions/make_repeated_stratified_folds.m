function folds = make_repeated_stratified_folds(strata,k,n_repeats,seed)
% Deterministic repeated stratified K-fold labels.
strata = string(strata(:));
n = numel(strata);
folds = zeros(n,n_repeats);
levels = unique(strata,'stable');
for r = 1:n_repeats
    rng(seed+r,'twister');
    labels = zeros(n,1);
    for j = 1:numel(levels)
        idx = find(strata==levels(j));
        idx = idx(randperm(numel(idx)));
        labels(idx) = mod((0:numel(idx)-1)',k)+1;
    end
    folds(:,r) = labels;
end
end
