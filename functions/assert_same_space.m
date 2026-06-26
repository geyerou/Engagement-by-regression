function [ok, message] = assert_same_space(info_a, info_b, tolerance)
if nargin < 3, tolerance = 1e-4; end
same_size = isequal(info_a.ImageSize(1:3), info_b.ImageSize(1:3));
same_pix = max(abs(double(info_a.PixelDimensions(1:3)) - ...
    double(info_b.PixelDimensions(1:3)))) < tolerance;
Ta = double(info_a.Transform.T);
Tb = double(info_b.Transform.T);
same_affine = max(abs(Ta(:)-Tb(:))) < tolerance;
ok = same_size && same_pix && same_affine;
message = sprintf('size=%d; voxel=%d; affine=%d', ...
    same_size, same_pix, same_affine);
end
