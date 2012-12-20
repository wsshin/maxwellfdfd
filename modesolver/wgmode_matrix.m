function [A, En_Ht, gEh_HvEn, gEv_HhEn, Hn_Et] = wgmode_matrix(omega, eps_face_cell, mu_edge_cell, s_factor_cell, grid3d, normal_axis, intercept, pml_kind)
% The operator A acts on H fields, not E fields.  Therefore, A x = lambda x
% solves for H fields of an eigenmode of the system. For the eigenvalue
% equation, see Sec II of G. Veronis and S. Fan, Journal of Lightwave
% Technology, vol. 25, no. 9, pp.2511--2521.

chkarg(istypesizeof(omega, 'real') && omega > 0, '"omega" should be positive.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(eps_face_cell, 'complexcell', [1, Axis.count], grid3d.N), ...
	'"eps_face_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(mu_edge_cell, 'complexcell', [1, Axis.count], grid3d.N), ...
	'"mu_edge_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count, GK.count], [1 0]), ...
	'"s_factor_cell" should be %d-by-%d cell array whose each element is row vector with real elements.', Axis.count, GK.count);

chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');

chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
g = GK.dual;
ind_n = ismembc2(intercept, grid3d.l{normal_axis, g});  % ind_n == 0 if intercept is not in grid3d.l{normal_axis, g}
assert(ind_n ~= 0, ...
	'%s grid in %s-axis of "grid3d" does not have "intercept" as an element.', g, normal_axis);
chkarg(intercept > grid3d.lpml(normal_axis, Sign.n) && intercept < grid3d.lpml(normal_axis, Sign.p), ...
	'"intercept" should be outside PML in the "normal_axis" direction.')

chkarg(istypesizeof(pml_kind, 'PML'), '"pml_kind" should be instance of PML.');

grid2d = Grid2d(grid3d, normal_axis);
h = grid2d.axis(Dir.h);
v = grid2d.axis(Dir.v);
n = grid2d.normal_axis;

dims = int([h, v, n]);
eps_face_cell{Axis.x} = permute(eps_face_cell{Axis.x}, dims);
eps_face_cell{Axis.y} = permute(eps_face_cell{Axis.y}, dims);
eps_face_cell{Axis.z} = permute(eps_face_cell{Axis.z}, dims);

mu_edge_cell{Axis.x} = permute(mu_edge_cell{Axis.x}, dims);
mu_edge_cell{Axis.y} = permute(mu_edge_cell{Axis.y}, dims);
mu_edge_cell{Axis.z} = permute(mu_edge_cell{Axis.z}, dims);

% Get eps and mu at the intercept.
eps_hh = eps_face_cell{h}(:,:,ind_n);
eps_vv = eps_face_cell{v}(:,:,ind_n);
eps_nn = eps_face_cell{n}(:,:,ind_n);

mu_hh = mu_edge_cell{h}(:,:,ind_n);
mu_vv = mu_edge_cell{v}(:,:,ind_n);
mu_nn = mu_edge_cell{n}(:,:,ind_n);

% Check the homogeneity of eps and mu at "intercept" in the "normal_axis" direction.
% assert(ind_n >= 2);
% eps_nn_prev = eps_face_cell{n}(:,:,ind_n-1);
% chkarg(all(all(eps_nn == eps_nn_prev)), ...
% 	'"eps_edge_cell" is not homogeneous in the "normal_axis" direction at "intercept".')

% mu_hh_prev = mu_edge_cell{h}(:,:,ind_n-1);
% mu_vv_prev = mu_edge_cell{v}(:,:,ind_n-1);
% chkarg(all(all(mu_hh == mu_hh_prev)) && all(all(mu_vv == mu_vv_prev)), ...
% 	'"mu_face_cell" is not homogeneous in the "normal_axis" direction at "intercept".')

% Derive various parameters (for 2D problems) from grid2d.
bc = grid2d.bc;
Nh = grid2d.N(Dir.h); Nv = grid2d.N(Dir.v);
N = Nh*Nv;
dh_prim = grid2d.dl{Dir.h, GK.prim};
dv_prim = grid2d.dl{Dir.v, GK.prim};
dh_dual = grid2d.dl{Dir.h, GK.dual};
dv_dual = grid2d.dl{Dir.v, GK.dual};

% Set up the diagonal matrices from edge lengths.
[dh_prim, dv_prim] = ndgrid(dh_prim, dv_prim);
[dh_dual, dv_dual] = ndgrid(dh_dual, dv_dual);
% dh_prim = repmat(dh_prim.', [1 Nv]);
% dh_dual = repmat(dh_dual.', [1 Nv]);
% dv_prim = repmat(dv_prim, [Nh 1]);
% dv_dual = repmat(dv_dual, [Nh 1]);

% Set up s-parameters.
if isempty(s_factor_cell)
	sh_prim = ones(1, Nh);
	sh_dual = ones(1, Nh);
	sv_prim = ones(1, Nv);
	sv_dual = ones(1, Nv);
else
	sh_prim = s_factor_cell{h, GK.prim};
	sh_dual = s_factor_cell{h, GK.dual};
	sv_prim = s_factor_cell{v, GK.prim};
	sv_dual = s_factor_cell{v, GK.dual};
end

[sh_prim, sv_prim] = ndgrid(sh_prim, sv_prim);
[sh_dual, sv_dual] = ndgrid(sh_dual, sv_dual);
% sh_prim = repmat(sh_prim.', [1, Nv]);
% sh_dual = repmat(sh_dual.', [1, Nv]);
% sv_prim = repmat(sv_prim, [Nh, 1]);
% sv_dual = repmat(sv_dual, [Nh, 1]);

% Create eps and mu matrices.  Note that each of eps_pp, eps_qq, mu_rr can be
% either a scalar or matrix. 
if pml_kind == PML.u
	eps_hh = eps_hh .* sv_dual ./ sh_prim;  % eps_hh is multiplied to Eh, which is at primary grid location in h and at dual grid location in v.
	eps_vv = eps_vv .* sh_dual ./ sv_prim;  % eps_vv is multiplied to Ev, which is at dual grid location in h and at primary grip location in v.
	eps_nn = eps_nn .* sh_dual .* sv_dual;  % eps_nn is multiplied to En, which is at dual grid location in both h and v.
	mu_hh = mu_hh .* sv_prim ./ sh_dual;  % mu_hh is multiplied to Hh, which is at dual grid location in h and at primary grid location in v.
	mu_vv = mu_vv .* sh_prim ./ sv_dual;  % mu_vv is multiplied to Hv, which is at primary grid location in h and at dual grid location in v.
	mu_nn = mu_nn .* sh_prim .* sv_prim;  % mu_nn is multiplied to Hn, which is at primary grid location in both h and v.
else
	assert(pml_kind == PML.sc);
	dh_prim = dh_prim .* sh_prim;
	dh_dual = dh_dual .* sh_dual;
	dv_prim = dv_prim .* sv_prim;
	dv_dual = dv_dual .* sv_dual;
end

% Create differential operators.
DhHv = generate_DpHq(Dir.h, bc, dh_dual);
DvHh = generate_DpHq(Dir.v, bc, dv_dual);

DhHh = generate_DpHp(Dir.h, bc, dh_prim);
DvHv = generate_DpHp(Dir.v, bc, dv_prim);

DhEn = generate_DpEr(Dir.h, bc, dh_prim);
DvEn = generate_DpEr(Dir.v, bc, dv_prim);

DhHn = generate_DpHr(Dir.h, bc, dh_dual);
DvHn = generate_DpHr(Dir.v, bc, dv_dual);

DhEv = generate_DpEq(Dir.h, bc, dh_prim); 
DvEh = generate_DpEq(Dir.v, bc, dv_prim); 

% Create the operator.
%
% See Sec. II of G. Veronis and S. Fan, Journal of Lightwave Technology, vol.
% 25, no. 9, pp.2511--2521.
eps_vv_hh = [eps_vv(:); eps_hh(:)];
EPS_vv_hh = spdiags(eps_vv_hh, 0, 2*N, 2*N);
INV_EPS_nn = spdiags(1./eps_nn(:), 0, N, N);  % when eps_nn has Inf, "EPS_nn \ Mat" complains about singularity
INV_MU_nn = spdiags(1./mu_nn(:), 0, N, N);
MU_hh_vv = spdiags([mu_hh(:); mu_vv(:)], 0, 2*N, 2*N);
MU_hh = spdiags(mu_hh(:), 0, N, N);
MU_vv = spdiags(mu_vv(:), 0, N, N);
MU_nn = spdiags(mu_nn(:), 0, N, N);

A = -omega^2 * EPS_vv_hh * MU_hh_vv ...
    + EPS_vv_hh * [DvEn; -DhEn] * (INV_EPS_nn * [-DvHh, DhHv]) ...
    - [DhHn; DvHn] * (MU_nn \ ([DhHh, DvHv] * MU_hh_vv));

% ez_array = (DxHy*hy.array(:) - DyHx*hx.array(:))./eps_zz(:)/omega/sqrt(-1);
% ex_array = (sqrt(-1)*omega*(mu_yy(:).*hy.array(:)) - DxEz*ez.array(:))./gamma;
% ey_array = (-sqrt(-1)*omega*(mu_xx(:).*hx.array(:)) - DyEz*ez.array(:))./gamma;

En_Ht = (1/(omega*1i)) * INV_EPS_nn * [-DvHh DhHv];  % En from Ht
Hn_Et = (-1/(omega*1i)) * INV_MU_nn * [-DvEh DhEv];  % Hn from Et
gEh_HvEn = [1i*omega*MU_vv -DhEn];  % gamma*Eh from Hv, En
gEv_HhEn = [-1i*omega*MU_hh -DvEn];  % gamma*Ev from Hh, En

% Z = sparse(N,N);
% I = speye(N);
% gEt_Ht = [gEh_HvEn Z Z; Z Z gEv_HhEn] * [Z I; En_Ht; I Z; En_Ht];


% We should make sure that eigenvectors' elements corresponding to infinitely
% large elements of eps_vv_hh are zero.  If the corresponding rows and columns
% are zeros except for the diagonal elements, and if the diagonal elements are
% not equal to the eigenvalue, then the relevant elements of the eigenvectors
% are zeros.  The eigenvalue is gamma^2 = (1/L + 1i*beta)^2 = (1/L)^2 - beta^2 +
% 2i*beta.  Should it be real, beta = 0 and the eigenvalue (1/L)^2 cannot be
% negative.
ind_pec = isinf(abs(eps_vv_hh));
pec_mask = ones(size(ind_pec));
pec_mask(ind_pec) = 0;
PM = spdiags(pec_mask, 0, length(pec_mask), length(pec_mask));

diagvec = full(diag(A));
diagvec(ind_pec) = -1;
A = PM * A * PM;
A = spdiags(diagvec, 0, A);

% % To force the H field components normal to the Et = 0 boundary to be 0, mask
% % the matrix appropriately.
% mask_hx = ones(Nh,Nv);
% mask_hy = ones(Nh,Nv);
% 
% if bc(Dir.h, Sign.n) == BC.e
%     mask_hx(1,:) = 0;
% end
% 
% if bc(Dir.v, Sign.n) == BC.e
%     mask_hy(:,1) = 0;
% end
% 
% Mask = spdiags([mask_hx(:); mask_hy(:)], 0, 2*N, 2*N);
% A = Mask*A;

% Reorder the indices of the elements of A to reduce the bandwidth of A.
r = reordering_indices(Dir.count, [Nh Nv]);

A = A(r,r);



