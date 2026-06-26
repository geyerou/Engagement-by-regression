function labels = read_schaefer_labels(cfg)
T = readtable(cfg.schaefer_label_file, 'TextType', 'string');
T = T(T.index >= 1 & T.index <= 400, :);
T = sortrows(T, 'index');
if height(T) ~= 400 || any(T.index ~= (1:400)')
    error('Schaefer label CSV must contain exactly indices 1:400.');
end

net17_names = unique(T.yeo17_network, 'stable');
[~, net17_id] = ismember(T.yeo17_network, net17_names);

labels = struct();
labels.table = T;
labels.roi_names = T.name;
labels.hemisphere = T.hemisphere;
labels.network17_names = net17_names;
labels.network17_id = net17_id;
end
