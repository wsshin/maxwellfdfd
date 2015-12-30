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

% Below, when eps is the primary material parameter and mu is the dual material
% parameter, they are averaged as Mf_prim * eps * Mb_prim and Mb_dual * mu *
% Mf_dual.

% Choice 1 of averaging for obtaining a symmetric operator: line averaging for
% averaging along the field directions; arithmetic averating for averaging along
% the directions normal to the field directions.
Mf_prim = Mf;  % forward arithmetic averaging for interpolated primary fields at primary grid cell vertices
Mb_dual = Mb;  % backward arithmetic averaging for interpolated dual fields at dual grid cell vertices (= primary grid cell centers)
Mb_prim = cell(1, Axis.count);  % backward line averaging for primary fields
Mf_dual = cell(1, Axis.count);  % forward line averaging for dual fields
for w = Axis.elems
	Mb_prim{w} = create_spdiag(dl{GT.prim, w}.^-1) * Mb{w} * create_spdiag(dl{GT.dual, w});
	Mf_dual{w} = create_spdiag(dl{GT.dual, w}.^-1) * Mf{w} * create_spdiag(dl{GT.prim, w});
end

% % Choice 2 of averaging for obtaining a symmetric operator: line averaging for
% % backward averaging; arithmetic averating for forward averaging.
% Mf_prim = Mf;  % arithmetic averaging for forward averaging
% Mb_prim = cell(1, Axis.count);  % line averaging for backward averaging
% for w = Axis.elems
% 	Mb_prim{w} = create_spdiag(dl{GT.prim, w}.^-1) * Mb{w} * create_spdiag(dl{GT.dual, w});
% end
% Mf_dual = Mf_prim;
% Mb_dual = Mb_prim;

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
	[Mf_prim{Axis.x}, Z, Z; ...
	Z, Mf_prim{Axis.y}, Z; ...
	Z, Z, Mf_prim{Axis.z}] * ...
		prim_mat_diag{2} * ...
		[Z, Mb_prim{Axis.y}, Z; ...
		Z, Z, Mb_prim{Axis.z}; ...
		Mb_prim{Axis.x}, Z, Z] + ...
	[Mf_prim{Axis.x}, Z, Z; ...
	Z, Mf_prim{Axis.y}, Z; ...
	Z, Z, Mf_prim{Axis.z}] * ...
		prim_mat_diag{3} * ...
		[Z, Z, Mb_prim{Axis.z}; ...
		Mb_prim{Axis.x}, Z, Z; ...
		Z, Mb_prim{Axis.y}, Z];

% Mask the matrix before adding the diagonal entries.  This ensures the masked
% fields to be zero in the solution, as long as the source there is zero.
PRIM_MAT(:, ind_mask_p) = 0;
PRIM_MAT(ind_mask_p, :) = 0;

PRIM_MAT = prim_mat_diag{1} + PRIM_MAT;

DUAL_MAT = ...
	[Mb_dual{Axis.x}, Z, Z; ...
	Z, Mb_dual{Axis.y}, Z; ...
	Z, Z, Mb_dual{Axis.z}] * ...
		dual_mat_diag{2} * ...
		[Z, Mf_dual{Axis.y}, Z; ...
		Z, Z, Mf_dual{Axis.z}; ...
		Mf_dual{Axis.x}, Z, Z] + ...
	[Mb_dual{Axis.x}, Z, Z; ...
	Z, Mb_dual{Axis.y}, Z; ...
	Z, Z, Mb_dual{Axis.z}] * ...
		dual_mat_diag{3} * ...
		[Z, Z, Mf_dual{Axis.z}; ...
		Mf_dual{Axis.x}, Z, Z; ...
		Z, Mf_dual{Axis.y}, Z];

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
