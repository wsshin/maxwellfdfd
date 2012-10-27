function phi = flux_patch(scalar2d, rect)

chkarg(istypesizeof(scalar2d, 'Scalar2d'), '"scalar2d" should be instance of Scalar2d.');

grid2d = scalar2d.grid2d;
if nargin < 2  % no rect
	rect = grid2d.bound_plot(false);  % false: do not include PML
end
chkarg(istypesizeof(rect, 'real', [2 2]), '"rect" should be %d-by-%d array with real elements.');

chkarg(issorted(rect(Dir.h,:)) && issorted(rect(Dir.v,:)), 'each row of "rect" should be sorted in ascending order.');
chkarg(grid2d.contains(rect.'), '"rect" should be contained in grid.');

l = grid2d.lall(:, GK.prim);
[Xh, Yv] = meshgrid(l{:});

li = cell(1, Dir.count);
for d = Dir.elems
	ind = NaN(1, Sign.count);
	for s = Sign.elems
		bound = rect(d, s);
		i = ismembc2(bound, l{d});
		if i == 0  % bound is not in l{d}
			warning('FDS:interp', ['patch boundary %s = %f is not aligned with primary grid;', ...
				'Poynting vector is interpolated at the boundary and flux will be slightly inaccurate.'], ...
				char(grid2d.axis(d)), bound);
		end
			
		if s == Sign.n
			ind(s) = find(l{d} > bound, 1, 'first');
		else  % s == Sign.p
			ind(s) = find(l{d} < bound, 1, 'last');
		end
	end
	li{d} = [rect(d,Sign.n), l{d}(ind(Sign.n):ind(Sign.p)), rect(d,Sign.p)];
end
[XIh, YIv] = meshgrid(li{:});

C = scalar2d.array;
C = permute(C, int([Dir.v, Dir.h]));
CI = interp2(Xh, Yv, C, XIh, YIv);
CI = ipermute(CI, int([Dir.v, Dir.h]));

phi = trapz(li{Dir.h}, CI, int(Dir.h));  % integrate along horizontal direction
phi = trapz(li{Dir.v}, phi);  % integrate along vertical direction
