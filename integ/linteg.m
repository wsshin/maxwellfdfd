function val = linteg(field, pt_start, pt_end)
chkarg(istypesizeof(field, 'Scalar2d') || istypesizeof(field, 'Scalar3d'), ...
	'"field" should be instance of Scalard2d or Scalar3d.');

if istypesizeof(field, 'Scalar3d')
	grid = field.grid3d;
	dir = Axis.x;  % temporary direction
else
	grid = field.grid2d;
	dir = Dir.h;  % temporary direction
end
chkarg(istypesizeof(pt_start, 'real', [1 dir.count]), ...
	'"pt_start" should be length-%d row vector with real elements.', dir.count);
chkarg(istypesizeof(pt_end, 'real', [1 dir.count]), ...
	'"pt_start" should be length-%d row vector with real elements.', dir.count);

ind_dir = find(pt_start ~= pt_end);
chkarg(length(ind_dir) <= 1, 'line connecting "pt_start" and "pt_end" should be along Cartesian direction.');
if isempty(ind_dir)  % all(pt_start == pt_end)
	val = 0;
else
	dir = dir.elems(ind_dir);  % actual direction

	gt = field.gt_array;
	li = cell(1, dir.count);
	for d = setdiff(dir.elems, dir)  % direction orthogonal to line of integral
		l = grid.lall{d, gt(d)};
		coord = pt_start(d);
		i = ismembc2(coord, l);
		if i == 0  % bound is not in lall
			warning('Maxwell:interp', ['integral limit %s = %s is not aligned with %s grid;', ...
				'"field" value is interpolated at the boundary and line integral will be slightly inaccurate.'], ...
				char(d), num2str(coord), char(gt(d)));
		end

		li{d} = coord;
	end
	

	r = grid.lall{dir, gt(dir)};
	ra = grid.lall{dir, alter(gt(dir))};
	ind_array = NaN(1, Sign.count);
	ind_dr = NaN(1, Sign.count);
	limit = [pt_start(dir) pt_end(dir)];
	limit = sort(limit);  % pt_start(dir) could be larger than pt_end(dir)
	for s = Sign.elems
		bound = limit(s);
		if s == Sign.n
			ind_array(s) = find(r >= bound, 1, 'first');
			ind_dr(s) = find(ra > bound, 1, 'first');
		else  % s == Sign.p
			ind_array(s) = find(r <= bound, 1, 'last');
			ind_dr(s) = find(ra < bound, 1, 'last');
		end
	end
	
	li{dir} = r(ind_array(Sign.n):ind_array(Sign.p));
	dr = [limit(Sign.n), ra(ind_dr(Sign.n):ind_dr(Sign.p)), limit(Sign.p)];
	if pt_start(dir) > pt_end(dir)
		li{dir} = fliplr(li{dir});
		dr = fliplr(dr);
	end
	dr = diff(dr);
	
	l = grid.lall(dir.elems + dir.count*subsindex(gt));
	if istypesizeof(field, 'Scalar3d')
		[X, Y, Z] = ndgrid(l{:});
		[XI, YI, ZI] = ndgrid(li{:});

		V = field.array;
		VI = interpn(X, Y, Z, V, XI, YI, ZI);
		
		vals = VI(:);
	else  % field is Scalar2d
		[Xh, Yv] = ndgrid(l{:});
		[XIh, YIv] = ndgrid(li{:});
		
		C = field.array;
		CI = interpn(Xh, Yv, C, XIh, YIv);
		
		vals = CI(:);
	end
	
	val = dr * vals;
end
	
