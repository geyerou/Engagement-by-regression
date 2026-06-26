function write_map_like(values, voxel_indices, reference_file, output_file)
info = niftiinfo(reference_file);
vol = zeros(info.ImageSize(1:3), 'single');
vol(voxel_indices) = single(values);
out_dir = fileparts(output_file);
ensure_dir(out_dir);

[~,name,ext] = fileparts(output_file);
if strcmpi(ext,'.gz')
    [~,name2,~] = fileparts(name);
    nii_base = fullfile(out_dir, name2);
    compressed = true;
else
    nii_base = fullfile(out_dir, name);
    compressed = false;
end
niftiwrite(vol, nii_base, info, 'Compressed', compressed);
end
