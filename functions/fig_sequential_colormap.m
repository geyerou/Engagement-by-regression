function cmap = fig_sequential_colormap(n,scheme)
if nargin < 1 || isempty(n), n = 256; end
if nargin < 2 || isempty(scheme), scheme = 'bluegreen'; end
switch lower(string(scheme))
    case "bluegreen"
        anchors = [ ...
            0.060 0.110 0.210
            0.075 0.270 0.420
            0.000 0.475 0.416
            0.460 0.690 0.360
            0.955 0.815 0.310];
    case "rose"
        anchors = [ ...
            0.150 0.090 0.160
            0.350 0.160 0.300
            0.640 0.210 0.330
            0.850 0.420 0.330
            0.970 0.760 0.520];
    otherwise
        anchors = parula(6);
end
x = linspace(0,1,size(anchors,1));
xi = linspace(0,1,n);
cmap = interp1(x,anchors,xi,'pchip');
cmap = max(0,min(1,cmap));
end
