function scalar2d = slice_scalar3d(scalar3d, normal_axis, intercept)

chkarg(istypesizeof(scalar3d, 'Scalar3d'), '"scalar3d" should be instance of Scalar3d.');
chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');

grid3d = scalar3d.grid3d;
[h, v, n] = cycle(normal_axis);
chkarg(grid3d.comp(n).contains(intercept), ...
	'"intercept" should be contained in %s-axis range.', char(n));

hall = grid3d.lall{h, GK.prim};
vall = grid3d.lall{v, GK.prim};
nall = grid3d.lall{n, GK.prim};

in = ismembc2(intercept, nall);
if in == 0
	warning('FDS:interp', 'slice at %s = %f is not primary grid plane; fields are interpolated.', char(n), intercept);
end

li = {hall, vall, intercept};
[XIh, YIv, ZIn] = meshgrid(li{:});

V = scalar3d.array;
V = permute(V, int([v h n]));  % [v h n] rather than [h v n]
l = {hall, vall, nall};
[Xh, Yv, Zn] = meshgrid(l{:});

array = interp3(Xh, Yv, Zn, V, XIh, YIv, ZIn);
array = ipermute(array, int([Dir.v Dir.h]));
assert(ndims(array) == Dir.count);

grid2d = Grid2d(grid3d, n);
scalar2d = Scalar2d(array, grid2d, scalar3d.osc, scalar3d.physQcell, scalar3d.name, intercept);
