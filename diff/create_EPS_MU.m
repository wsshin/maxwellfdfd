function [EPS, MU] = create_EPS_MU(eps_array, mu_array, ge, dl_factor_cell, grid3d)

chkarg(istypesizeof(eps_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
	'"eps_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(mu_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
	'"mu_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for the E-field
chkarg(isempty(dl_factor_cell) || istypesizeof(dl_factor_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
	'"dl_factor_cell" should be empty, or %d-by-%d cell array whose each element is row vector with real elements.', Axis.count, GT.count);
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

eps_diag = create_mat_diag(eps_array, grid3d);
mu_diag = create_mat_diag(mu_array, grid3d);

%% Create Df and Db.
Mf = create_Ms(Sign.p, ge, grid3d);
Mb = create_Ms(Sign.n, ge, grid3d);

dl = cell(GT.count, Axis.count);  % defining dl = cell(Axis.count, GT.count) and using [dl{:,g}] below doesn't work
for g = GT.elems
	if isempty(dl_factor_cell)
		[dl{g,:}] = ndgrid(grid3d.dl{:,g});
	else
		dl_cell = mult_vec(dl_factor_cell(:,g), grid3d.dl(:,g));
		[dl{g,:}] = ndgrid(dl_cell{:});
	end
end

Mfp = Mf;  % forward averaging for primary fields; it is simply arithmetic averaging
Mbd = Mb;  % backward averaging for dual fields; it is simply arithmetic averaging

Mbp = cell(1, Axis.count);  % backward averaging for primary fields; it is line averaging
Mfd = cell(1, Axis.count);  % forward averaging for dual fields; it is line averaging
for w = Axis.elems
	Mbp{w} = create_spdiag(dl{GT.dual, w}.^-1) * Mb{w} * create_spdiag(dl{GT.prim, w});
	Mfd{w} = create_spdiag(dl{GT.prim, w}.^-1) * Mf{w} * create_spdiag(dl{GT.dual, w});
end

%% Create mask matrices.
[ind_mask_p, ind_mask_d] = create_masks(ge, grid3d);

%% Form averaging matrices.
Ntot = prod(grid3d.N);
Z = sparse(Ntot, Ntot);

if ge == GT.prim
	prim_mat_diag = eps_diag;
	dual_mat_diag = mu_diag;
else  % ge == GT.dual
	prim_mat_diag = mu_diag;
	dual_mat_diag = eps_diag;
end

PRIM_MAT = ...
	[Mfp{Axis.x}, Z, Z; ...
	Z, Mfp{Axis.y}, Z; ...
	Z, Z, Mfp{Axis.z}] * ...
		prim_mat_diag{2} * ...
		[Z, Mbp{Axis.y}, Z; ...
		Z, Z, Mbp{Axis.z}; ...
		Mbp{Axis.x}, Z, Z] + ...
	[Mfp{Axis.x}, Z, Z; ...
	Z, Mfp{Axis.y}, Z; ...
	Z, Z, Mfp{Axis.z}] * ...
		prim_mat_diag{3} * ...
		[Z, Z, Mbp{Axis.z}; ...
		Mbp{Axis.x}, Z, Z; ...
		Z, Mbp{Axis.y}, Z];

% Mask the matrix before adding the diagonal entries.  This ensures the masked
% fields to be zero in the solution, as long as the source there is zero.
PRIM_MAT(:, ind_mask_p) = 0;
PRIM_MAT(ind_mask_p, :) = 0;

PRIM_MAT = prim_mat_diag{1} + PRIM_MAT;

DUAL_MAT = ...
	[Mbd{Axis.x}, Z, Z; ...
	Z, Mbd{Axis.y}, Z; ...
	Z, Z, Mbd{Axis.z}] * ...
		dual_mat_diag{2} * ...
		[Z, Mfd{Axis.y}, Z; ...
		Z, Z, Mfd{Axis.z}; ...
		Mfd{Axis.x}, Z, Z] + ...
	[Mbd{Axis.x}, Z, Z; ...
	Z, Mbd{Axis.y}, Z; ...
	Z, Z, Mbd{Axis.z}] * ...
		dual_mat_diag{3} * ...
		[Z, Z, Mfd{Axis.z}; ...
		Mfd{Axis.x}, Z, Z; ...
		Z, Mfd{Axis.y}, Z];

% Mask the matrix before adding the diagonal entries.  This ensures the masked
% fields to be zero in the solution, as long as the source there is zero.
DUAL_MAT(:, ind_mask_d) = 0;
DUAL_MAT(ind_mask_d, :) = 0;

DUAL_MAT = dual_mat_diag{1} + DUAL_MAT;

if ge == GT.prim
	EPS = PRIM_MAT;
	MU = DUAL_MAT;
else  % ge == GT.dual
	EPS = DUAL_MAT;
	MU = PRIM_MAT;
end
