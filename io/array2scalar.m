function scalar = array2scalar(F_array, physQ, grid, axis, gk, osc, intercept)
chkarg(istypesizeof(physQ, 'PhysQ'), '"physQ" should be instance of PhysQ.');
chkarg(istypesizeof(gk, 'GK'), '"gk" should be instance of GK.');  % GK.dual for E-like field; GK.prim for H-like field
chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instancef of Oscillation.');

chkarg(istypesizeof(grid, 'Grid2d') || istypesizeof(grid, 'Grid3d'), ...
	'"grid" should be instance of Grid2d or Grid3d.');
is3D = true;
if istypesizeof(grid, 'Grid2d')
	is3D = false;
end

if is3D
	grid3d = grid;
	chkarg(istypesizeof(F_array, 'complex', grid3d.N), ...
		'"F_array" should be %d-by-%d-by-%d array with complex elements.', ...
		grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));

	l = grid3d.lall(:, gk);
	l{axis} = grid3d.lall{axis, alter(gk)};
	[X, Y, Z] = meshgrid(l{:});

	V = F_array;
	for w = Axis.elems
		V = attach_extra_F(V, gk, axis, grid3d.bc(w,Sign.n), w);
	end
	V = permute(V, int([Axis.y, Axis.x, Axis.z]));
	
	li = grid3d.lall(:,GK.prim);
	[XI, YI, ZI] = meshgrid(li{:});
	
	VI = interp3(X, Y, Z, V, XI, YI, ZI);
	VI = ipermute(VI, int([Axis.y, Axis.x, Axis.z]));
	scalar = Scalar3d(VI, grid3d, osc, physQ, [physQ.symbol, '_', char(axis)]);
else  % grid is Grid2d
	if nargin < 7  % no intercept
		intercept = NaN;
	end
	grid2d = grid;
	chkarg(istypesizeof(F_array, 'complex', grid2d.N), ...
		'"F_array" should be %d-by-%d array with complex elements.', grid2d.N(Dir.h), grid2d.N(Dir.v));

	l = grid2d.lall(:, gk);
	axis_temp = Axis.z;
	if axis == grid2d.axis(Dir.h)
		l{Dir.h} = grid2d.lall{Dir.h, alter(gk)};
		axis_temp = Axis.x;
	elseif axis == grid2d.axis(Dir.v)
		l{Dir.v} = grid2d.lall{Dir.v, alter(gk)};
		axis_temp = Axis.y;
	end	
	[Xh, Yv] = meshgrid(l{:});

	C = F_array;
	C = attach_extra_F(C, gk, axis_temp, grid2d.bc(Dir.h,Sign.n), Axis.x);  % treat C as 3D array with one element in 3rd dimension
	C = attach_extra_F(C, gk, axis_temp, grid2d.bc(Dir.v,Sign.n), Axis.y);  % treat C as 3D array with one element in 3rd dimension
	C = permute(C, int([Dir.v, Dir.h]));

	li = grid2d.lall(:,GK.prim);
	[XIh, YIv] = meshgrid(li{:});

	CI = interp2(Xh, Yv, C, XIh, YIv);
	CI = ipermute(CI, int([Dir.v, Dir.h]));
	scalar = Scalar2d(CI, grid2d, osc, physQ, [physQ.symbol, '_', char(axis)], intercept);
end


function Fw_array = attach_extra_F(Fw_array, gk, w, bc_vn, v)
% Augment the Fw array with ghost boundary elements in the v-axis.
% Surprisingly, this function combines attach_extra_E() and attach_extra_H().

chkarg(istypesizeof(Fw_array, 'complex', [0 0 0]), '"Fw_array" should be 3D array with complex elements.');
chkarg(istypesizeof(gk, 'GK'), '"physQ" should be instance of GK.');
chkarg(istypesizeof(w, 'Axis'), '"w" should be instance of Axis.');
chkarg(istypesizeof(bc_vn, 'BC'), '"bc_vn" should be instance of BC.');
chkarg(istypesizeof(v, 'Axis'), '"v" should be instance of Axis.');

ind_n = {':',':',':'};
ind_p = {':',':',':'};
Nv = size(Fw_array, int(v));
if (gk == GK.dual && v == w) || (gk == GK.prim && v ~= w)  % E: gk == GK.dual, H: gk == GK.prim
	if bc_vn == BC.p
		ind_p{v} = 1;
		Fw_array = cat(int(v), Fw_array, Fw_array(ind_p{:}));
	else  % bc_vp == BC.m
		size_extra = size(Fw_array);
		size_extra(v) = 1;
		Fw_array = cat(int(v), Fw_array, zeros(size_extra));
	end
else  % (gk == GK.dual && v ~= w) || (gk == GK.prim && v == w)
	if bc_vn == BC.p
		ind_n{v} = Nv;
		ind_p{v} = 1;
		Fw_array = cat(int(v), Fw_array(ind_n{:}), Fw_array, Fw_array(ind_p{:}));
	else  % bc_vn == BC.m || BC.e, bc_vp == BC.m
		ind_n{v} = 1;
		ind_p{v} = Nv;
		if bc_vn == BC.e
			Fw_array = cat(int(v), -Fw_array(ind_n{:}), Fw_array, Fw_array(ind_p{:}));
		else  % bc_vn == BC.m
			Fw_array = cat(int(v), Fw_array(ind_n{:}), Fw_array, Fw_array(ind_p{:}));
		end
	end
end
