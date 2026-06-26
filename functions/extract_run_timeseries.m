function run_data = extract_run_timeseries(func_file, atlas_file, wm_file, cfg)
func_info = niftiinfo(func_file);
atlas_info = niftiinfo(atlas_file);
wm_info = niftiinfo(wm_file);
[ok1,msg1] = assert_same_space(func_info, atlas_info);
[ok2,msg2] = assert_same_space(func_info, wm_info);
if ~ok1 || ~ok2
    error('Space mismatch: atlas %s; WM %s', msg1, msg2);
end

func = single(niftiread(func_info));
atlas = round(single(niftiread(atlas_info)));
wm = niftiread(wm_info) > 0;
T = size(func,4);
wm_idx = find(wm);

X = zeros(T,400,'single');
flat = reshape(func, [], T);
for roi = 1:400
    idx = find(atlas == roi);
    if isempty(idx)
        error('ROI %d has no voxels.', roi);
    end
    X(:,roi) = mean(flat(idx,:), 1, 'omitnan')';
end

if cfg.do_wm_mask_smoothing
    if ~license('test','image_toolbox')
        error('Mask-normalized smoothing requires Image Processing Toolbox.');
    end
    sigma_vox = cfg.wm_smooth_fwhm_mm / ...
        (2*sqrt(2*log(2))*mean(func_info.PixelDimensions(1:3)));
    den = imgaussfilt3(single(wm), sigma_vox, 'Padding', 0);
    Y = zeros(T,numel(wm_idx),'single');
    for t = 1:T
        vol = func(:,:,:,t);
        num = imgaussfilt3(vol .* single(wm), sigma_vox, 'Padding', 0);
        smoothed = num ./ max(den, eps('single'));
        Y(t,:) = smoothed(wm_idx);
    end
else
    Y = flat(wm_idx,:)';
end
clear func flat

% The upstream preprocessing added each voxel's baseline back after filtering.
% Remove run means here so the model targets BOLD fluctuations, not baselines.
X = X - mean(X,1,'omitnan');
Y = Y - mean(Y,1,'omitnan');
X(~isfinite(X)) = 0;
Y(~isfinite(Y)) = 0;

run_data = struct();
run_data.X_gm_raw = X;
run_data.Y_wm_raw = Y;
run_data.wm_voxel_indices = wm_idx;
run_data.TR = double(func_info.PixelDimensions(4));
run_data.n_timepoints = T;
run_data.source_file = func_file;
run_data.wm_smoothed = cfg.do_wm_mask_smoothing;
end
