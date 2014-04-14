function [A, En_Htr, gEh_HvEn, gEv_HhEn, Hn_Etr] = wgmode_matrix(ge, pml, omega, eps_cell, mu_cell, s_factor_cell, grid3d, normal_axis, intercept)
% The operator A acts on H fields, not E fields.  Therefore, A x = lambda x
% solves for H fields of an eigenmode of the system. For the eigenvalue
% equation, see Sec II of G. Veronis and S. Fan, Journal of Lightwave
% Technology, vol. 25, no. 9, pp.2511--2521.

chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');
chkarg(istypesizeof(pml, 'PML'), '"pml" should be instance of PML.');
chkarg(istypesizeof(omega, 'real') && omega > 0, '"omega" should be positive.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(eps_cell, 'complexcell', [1, Axis.count], grid3d.N), ...
	'"eps_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(mu_cell, 'complexcell', [1, Axis.count], grid3d.N), ...
	'"mu_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count, GT.count], [1 0]), ...
	'"s_factor_cell" should be %d-by-%d cell array whose each element is row vector with real elements.', Axis.count, GT.count);

chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');

chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
ind_n = ismembc2(intercept, grid3d.l{normal_axis, ge});  % ind_n == 0 if intercept is not in grid3d.l{normal_axis, ge}
if ind_n == 0
	ind_n = ismembc2(intercept, grid3d.l{normal_axis, alter(ge)});
end
assert(ind_n ~= 0, ...
	'%s grid in %s-axis of "grid3d" does not have "intercept" as an element.', char(ge), char(normal_axis));
chkarg(intercept > grid3d.lpml(normal_axis, Sign.n) && intercept < grid3d.lpml(normal_axis, Sign.p), ...
	'"intercept" should be outside PML in the "normal_axis" direction.')


grid2d = Grid2d(grid3d, normal_axis);
h = grid2d.axis(Dir.h);
v = grid2d.axis(Dir.v);
n = grid2d.normal_axis;

g3d = Grid3d(grid3d.unit, ...
	{grid2d.lall{Dir.h,GT.prim}, grid2d.lall{Dir.v,GT.prim}, [0 1]}, ...
	[grid2d.Npml; 0 0], [grid2d.bc, BC.p]);

dims = int([h, v, n]);
eps_cell{Axis.x} = permute(eps_cell{Axis.x}, dims);
eps_cell{Axis.y} = permute(eps_cell{Axis.y}, dims);
eps_cell{Axis.z} = permute(eps_cell{Axis.z}, dims);

mu_cell{Axis.x} = permute(mu_cell{Axis.x}, dims);
mu_cell{Axis.y} = permute(mu_cell{Axis.y}, dims);
mu_cell{Axis.z} = permute(mu_cell{Axis.z}, dims);

% Get eps and mu at the intercept.
eps_h = eps_cell{h}(:,:,ind_n);
eps_v = eps_cell{v}(:,:,ind_n);
eps_n = eps_cell{n}(:,:,ind_n);

mu_h = mu_cell{h}(:,:,ind_n);
mu_v = mu_cell{v}(:,:,ind_n);
mu_n = mu_cell{n}(:,:,ind_n);

% Check the homogeneity of eps and mu at "intercept" in the "normal_axis" direction.
% assert(ind_n >= 2);
% eps_nn_prev = eps_face_cell{n}(:,:,ind_n-1);
% chkarg(all(all(eps_nn == eps_nn_prev)), ...
% 	'"eps_edge_cell" is not homogeneous in the "normal_axis" direction at "intercept".')

% mu_hh_prev = mu_edge_cell{h}(:,:,ind_n-1);
% mu_vv_prev = mu_edge_cell{v}(:,:,ind_n-1);
% chkarg(all(all(mu_hh == mu_hh_prev)) && all(all(mu_vv == mu_vv_prev)), ...
% 	'"mu_face_cell" is not homogeneous in the "normal_axis" direction at "intercept".')

% Create dl_factor_cell.
dl_factor_cell = [];
if pml == PML.sc
	dl_factor_cell = s_factor_cell([h v],:);
end

% Create eps and mu matrices.
if pml == PML.u
	[smh, smv] = ndgrid(s_factor_cell{h, alter(ge)}, s_factor_cell{v, alter(ge)});
	[seh, sev] = ndgrid(s_factor_cell{h, ge}, s_factor_cell{v, ge});
	
	eps_h = eps_h .* sev ./ smh;
	eps_v = eps_v .* seh ./ smv;
	eps_n = eps_n .* seh .* sev;
	
	mu_h = mu_h .* smv ./ seh;
	mu_v = mu_v .* smh ./ sev;
	mu_n = mu_n .* smh .* smv;
end

% Create differential operators
if ge == GT.prim
	s = Sign.n;
else
	s = Sign.p;
end
Ds = create_Ds(s, ge, dl_factor_cell, grid2d);
Da = create_Ds(alter(s), ge, dl_factor_cell, grid2d);


% DhHh = create_Ds(s, ge, dl_cell, grid2d);
% DvHv = create_Ds(s, ge, dl_cell, grid2d);
% 
% DhHv = create_Ds(alter(s), ge, dl_cell, grid2d);
% DvHh = create_Ds(alter(s), ge, dl_cell, grid2d);
% 
% DhEn = create_Ds(s, ge, dl_cell, grid2d);
% DvEn = create_Ds(s, ge, dl_cell, grid2d);
% 
% DhHn = create_Ds(alter(s), ge, dl_cell, grid2d);
% DvHn = create_Ds(alter(s), ge, dl_cell, grid2d);
% 
% DhEv = create_Ds(s, ge, dl_cell, grid2d);
% DvEh = create_Ds(s, ge, dl_cell, grid2d);



% Create the operator.
%
% See Sec. II of G. Veronis and S. Fan, Journal of Lightwave Technology, vol.
% 25, no. 9, pp.2511--2521.
eps_vh = [eps_v(:); eps_h(:)];
EPS_vh = create_spdiag(eps_vh);
INV_EPS_n = create_spdiag(1./eps_n(:));  % avoid ill-defined "EPS_nn \ Mat" operation when eps_nn has Inf
INV_MU_n = create_spdiag(1./mu_n(:));
MU_hv = create_spdiag([mu_h(:); mu_v(:)]);
MU_h = create_spdiag(mu_h(:));
MU_v = create_spdiag(mu_v(:));
MU_n = create_spdiag(mu_n(:));

if ge == GT.prim
	[ind_ME, ind_MH] = create_masks(ge, g3d);
else
	[ind_MH, ind_ME] = create_masks(ge, g3d);
end

ind_MH = reshape(ind_MH, [], Axis.count);
ind_MHh = ind_MH(:, Axis.x);
ind_MHv = ind_MH(:, Axis.y);
ind_MHn = ind_MH(:, Axis.z);

ind_ME = reshape(ind_ME, [], Axis.count);
ind_MEh = ind_ME(:, Axis.x);
ind_MEv = ind_ME(:, Axis.y);
ind_MEn = ind_ME(:, Axis.z);

DhHh = Da{Dir.h}; DhHh(:,ind_MHh) = 0;
DvHv = Da{Dir.v}; DvHv(:,ind_MHv) = 0;

DhHv = Ds{Dir.h}; DhHv(:,ind_MHv) = 0;
DvHh = Ds{Dir.v}; DvHh(:,ind_MHh) = 0;

DhEn = Da{Dir.h}; DhEn(:,ind_MEn) = 0;
DvEn = Da{Dir.v}; DvEn(:,ind_MEn) = 0;

DhHn = Ds{Dir.h}; DhHn(:,ind_MHn) = 0;
DvHn = Ds{Dir.v}; DvHn(:,ind_MHn) = 0;

DhEv = Da{Dir.h}; DhEv(:,ind_MEv) = 0;
DvEh = Da{Dir.v}; DvEh(:,ind_MEh) = 0;

A = -omega^2 * EPS_vh * MU_hv ...
    + EPS_vh * [DvEn; -DhEn] * (INV_EPS_n * [-DvHh, DhHv]) ...
    - [DhHn; DvHn] * (MU_n \ ([DhHh, DvHv] * MU_hv));
% A = -omega^2 * EPS_vh * MU_hv ...
%     + EPS_vh * [Da{Dir.v}; -Da{Dir.h}] * (INV_EPS_n * [-Ds{Dir.v}, Ds{Dir.h}]) ...
%     - [Ds{Dir.h}; Ds{Dir.v}] * (MU_n \ ([Da{Dir.h}, Da{Dir.v}] * MU_hv));

% ez_array = (DxHy*hy.array(:) - DyHx*hx.array(:))./eps_zz(:)/omega/sqrt(-1);
% ex_array = (sqrt(-1)*omega*(mu_yy(:).*hy.array(:)) - DxEz*ez.array(:))./gamma;
% ey_array = (-sqrt(-1)*omega*(mu_xx(:).*hx.array(:)) - DyEz*ez.array(:))./gamma;

En_Htr = (1/(omega*1i)) * INV_EPS_n * [-DvHh DhHv]; En_Htr(ind_MEn,:) = 0;  % En from Htr
Hn_Etr = (-1/(omega*1i)) * INV_MU_n * [-DvEh DhEv]; Hn_Etr(ind_MHn,:) = 0;  % Hn from Etr
gEh_HvEn = [1i*omega*MU_v -DhEn]; gEh_HvEn(ind_MEh,:) = 0;  % gamma*Eh from Hv, En
gEv_HhEn = [-1i*omega*MU_h -DvEn]; gEv_HhEn(ind_MEv,:) = 0;  % gamma*Ev from Hh, En
% En_Htr = (1/(omega*1i)) * INV_EPS_n * [-Ds{Dir.v}, Ds{Dir.h}];  % En from Htr
% Hn_Etr = (-1/(omega*1i)) * INV_MU_n * [-Da{Dir.v}, Da{Dir.h}];  % Hn from Etr
% gEh_HvEn = [1i*omega*MU_v -Da{Dir.h}];  % gamma*Eh from Hv, En
% gEv_HhEn = [-1i*omega*MU_h -Da{Dir.v}];  % gamma*Ev from Hh, En


% Z = sparse(N,N);
% I = speye(N);
% gEtr_Htr = [gEh_HvEn Z Z; Z Z gEv_HhEn] * [Z I; En_Htr; I Z; En_Htr];



% We should make sure that eigenvectors' elements corresponding to infinitely
% large elements of eps_vv_hh are zero.  If the corresponding rows and columns
% are zeros except for the diagonal elements, and if the corresponding diagonal
% elements are not equal to the eigenvalue, then the relevant elements of the
% eigenvectors are zeros.  The eigenvalue is gamma^2 = (1/L + 1i*beta)^2 =
% (1/L)^2 - beta^2 + 2i*beta.  Should it be real, beta = 0 and the eigenvalue
% (1/L)^2 cannot be negative, so we set the corresponding diagonal elements
% negative.
ind_pec = isinf(abs(eps_vh));
% pec_mask(ind_pec) = 0;
% PM = create_spdiag(pec_mask);

diagvec = full(diag(A));
diagvec(ind_pec) = -1;
% A = PM * A * PM;
A(:,ind_pec) = 0;
A(ind_pec,:) = 0;
A = spdiags(diagvec, 0, A);

% % To force the H field components normal to the Etr = 0 boundary to be 0, mask
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

ind_mask_h = [ind_MHh(:); ind_MHv(:)];
% A(ind_maskh, :) = 0;
% diagvec = full(diag(A));
% diagvec(ind_maskh) = -1;
% % A = PM * A * PM;
A(:,ind_mask_h) = 0;
A(ind_mask_h,:) = 0;
% A = spdiags(diagvec, 0, A);

% Reorder the indices of the elements of A to reduce the bandwidth of A.
r = reordering_indices(Dir.count, grid2d.N);

A = A(r,r);



