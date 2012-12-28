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

	li = cell(1, dir.count);
	for d = setdiff(dir.elems, dir)
		lall = grid.lall{d, GK.prim};
		bound = pt_start(d);
		i = ismembc2(bound, lall);
		if i == 0  % bound is not in lall
			warning('FDS:interp', ['integral limit %s = %f is not aligned with primary grid;', ...
				'"field" value is interpolated at the boundary and line integral will be slightly inaccurate.'], ...
				char(d), bound);
		end

		li{d} = bound;
	end
	

	lall = grid.lall{dir, GK.prim};
	ind = NaN(1, Sign.count);
	limit = [pt_start(dir) pt_end(dir)];
	limit = sort(limit);  % pt_start(dir) could be larger than pt_end(dir)
	for s = Sign.elems
		bound = limit(s);
		i = ismembc2(bound, lall);
		if i == 0  % bound is not in lall
			warning('FDS:interp', ['integral limit %s = %f is not aligned with primary grid;', ...
				'"field" value is interpolated at the boundary and line integral will be slightly inaccurate.'], ...
				char(dir), bound);
		end

		if s == Sign.n
			ind(s) = find(lall > bound, 1, 'first');
		else  % s == Sign.p
			ind(s) = find(lall < bound, 1, 'last');
		end
	end
	
	li{dir} = [limit(Sign.n), lall(ind(Sign.n):ind(Sign.p)), limit(Sign.p)];
	if pt_start(dir) > pt_end(dir)
		li{dir} = fliplr(li{dir});
	end
	
	l = grid.lall(:, GK.prim);
	if istypesizeof(field, 'Scalar3d')
		[X, Y, Z] = meshgrid(l{:});
		[XI, YI, ZI] = meshgrid(li{:});

		V = field.array;
		V = permute(V, int([Axis.y Axis.x Axis.z]));
		VI = interp3(X, Y, Z, V, XI, YI, ZI);
		VI = ipermute(VI, int([Axis.y Axis.x Axis.z]));
		
		vals = VI(:);
	else  % field is Scalar2d
		[Xh, Yv] = meshgrid(l{:});
		[XIh, YIv] = meshgrid(li{:});
		
		C = field.array;
		C = permute(C, int([Dir.v Dir.h]));
		CI = interp2(Xh, Yv, C, XIh, YIv);
		CI = ipermute(CI, int([Dir.v Dir.h]));
		
		vals = CI(:);
	end
	
	vals = vals.';  % make vals a row vector
	val = trapz(li{dir}, vals);
end
	
