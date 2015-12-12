function scalar = array2scalar(F_array, grid, ge, ft, axis, osc, physQ, intercept)

chkarg(istypesizeof(grid, 'Grid2d') || istypesizeof(grid, 'Grid3d'), ...
	'"grid" should be instance of Grid2d or Grid3d.');
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');
chkarg(istypesizeof(ft, 'FT'), '"ft" should be instance of FT.');
chkarg(istypesizeof(axis, 'Axis'), '"axis" should be instance of Axis.');
chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
chkarg(istypesizeof(physQ, 'PhysQ'), '"physQ" should be instance of PhysQ.');

is3D = true;
if istypesizeof(grid, 'Grid2d')
	is3D = false;
end

if ft == FT.e
    gt = ge;  % grid type for E-field
else  % ft == FT.h
    gt = alter(ge);  % grid type for H-field
end

gt_array = gt(ones(1, Axis.count));
gt_array(axis) = alter(gt);

if is3D
	grid3d = grid;
	chkarg(istypesizeof(F_array, 'arbitrary', grid3d.N), '"F_array" should be %d-by-%d-by-%d array.', grid3d.Ncell{:});

	V = F_array;
	for w = Axis.elems
		V = attach_extra_F(V, ft, gt_array(w), axis, grid3d.bc(w), w, grid3d.kBloch(w)*grid3d.L(w));
	end
	
	scalar = Scalar3d(V, grid3d, gt_array, osc, physQ, [physQ.symbol, '_', char(axis)]);
else  % grid is Grid2d
	if nargin < 7  % no intercept
		intercept = NaN;
	end
	chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');  % NaN is real
	
	grid2d = grid;
	chkarg(istypesizeof(F_array, 'arbitrary', grid2d.N), '"F_array" should be %d-by-%d array.', grid2d.Ncell{:});

	gt_array = gt_array(grid2d.axis);

	if axis == grid2d.axis(Dir.h)
		axis_temp = Axis.x;
	elseif axis == grid2d.axis(Dir.v)
		axis_temp = Axis.y;
	else  % axis == grid2d.normal_axis
		axis_temp = Axis.z;	
	end
	
	C = F_array;
	C = attach_extra_F(C, ft, gt_array(Dir.h), axis_temp, grid2d.bc(Dir.h), Axis.x, grid2d.kBloch(Dir.h)*grid2d.L(Dir.h));  % treat C as 3D array with one element in 3rd dimension
	C = attach_extra_F(C, ft, gt_array(Dir.v), axis_temp, grid2d.bc(Dir.v), Axis.y, grid2d.kBloch(Dir.v)*grid2d.L(Dir.v));  % treat C as 3D array with one element in 3rd dimension

	scalar = Scalar2d(C, grid2d, gt_array, osc, physQ, [physQ.symbol, '_', char(axis)], intercept);
end

function Fw_array = attach_extra_F(Fw_array, ft, gt, w, bc_v, v, kvLv)
% Augment the Fw array with ghost boundary elements in the v-axis.
% Surprisingly, this function combines attach_extra_E() and attach_extra_H().

% chkarg(istypesizeof(Fw_array, 'complex', [0 0 0]), '"Fw_array" should be 3D array with complex elements.');
chkarg(istypesizeof(ft, 'FT'), '"ft" should be instance of FT.');
chkarg(istypesizeof(gt, 'GT'), '"gt" should be instance of GT.');
chkarg(istypesizeof(w, 'Axis'), '"w" should be instance of Axis.');
chkarg(istypesizeof(bc_v, 'BC'), '"bc_v" should be instance of BC.');
chkarg(istypesizeof(v, 'Axis'), '"v" should be instance of Axis.');
chkarg(istypesizeof(kvLv, 'real'), '"kvLv" should be real.');

ind_n = {':',':',':'};
ind_p = {':',':',':'};
Nv = size(Fw_array, int(v));
if gt == GT.prim
	if bc_v == BC.p
		ind_p{v} = 1;
		Fw_array = cat(int(v), Fw_array, Fw_array(ind_p{:})*exp(-1i*kvLv));
	else  % bc_v == BC.e or BC.m
		% At the ghost point tangential primary fields are always zero, so
		% normal dual fields are also zero.
		size_extra = size(Fw_array);
		size_extra(v) = 1;
		Fw_array = cat(int(v), Fw_array, zeros(size_extra));
	end
else  % gt == GT.dual
	if bc_v == BC.p
		ind_n{v} = Nv;
		ind_p{v} = 1;
		Fw_array = cat(int(v), Fw_array(ind_n{:})*exp(1i*kvLv), Fw_array, Fw_array(ind_p{:})*exp(-1i*kvLv));
	else  % bc_v == BC.e or BC.m
		ind_n{v} = 1;
		ind_p{v} = Nv;
		ft_match_bc = (ft == FT.e && bc_v == BC.e) || (ft == FT.h && bc_v == BC.m);
		if (v ~= w && ft_match_bc) || (v == w && ~ft_match_bc)
			Fw_array = cat(int(v), -Fw_array(ind_n{:}), Fw_array, Fw_array(ind_p{:}));
		else  % (v == w && ft_match_bc) || (v ~= w && ~ft_match_bc)
			Fw_array = cat(int(v), Fw_array(ind_n{:}), Fw_array, Fw_array(ind_p{:}));
		end
	end
end
