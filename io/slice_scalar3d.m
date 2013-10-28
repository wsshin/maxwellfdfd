function scalar2d = slice_scalar3d(scalar3d, normal_axis, intercept)

chkarg(istypesizeof(scalar3d, 'Scalar3d'), '"scalar3d" should be instance of Scalar3d.');
chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');

grid3d = scalar3d.grid3d;
[h, v, n] = cycle(normal_axis);
chkarg(grid3d.comp(n).contains(intercept), ...
	'"intercept" should be contained in %s-axis range.', char(n));

lall = grid3d.lall(Axis.elems + Axis.count*subsindex(scalar3d.gt_array));
nall = lall{n};

ind = {':', ':', ':'};
indn = ismembc2(intercept, nall);
if indn ~= 0
	ind{n} = indn;
	array = scalar3d.array(ind{:});
else
% 	warning('Maxwell:interp', 'slice at %s = %s is not grid plane where %s is defined; fields are interpolated.', ...
% 		char(n), num2str(intercept), scalar3d.name);

	indn = find(nall < intercept, 1, 'last');
	ind{n} = [indn indn+1];

	V = scalar3d.array(ind{:});
	l = lall;
	l{n} = nall(ind{n});
	[Xh, Yv, Zn] = ndgrid(l{:});

	li = lall;
	li{n} = intercept;
	[XIh, YIv, ZIn] = ndgrid(li{:});

	array = interpn(Xh, Yv, Zn, V, XIh, YIv, ZIn);
end

array = permute(array, int([h v n]));
assert(ndims(array) == Dir.count);

grid2d = Grid2d(grid3d, n);
scalar2d = Scalar2d(array, grid2d, scalar3d.gt_array([h v]), scalar3d.osc, scalar3d.physQcell, scalar3d.name, intercept);
