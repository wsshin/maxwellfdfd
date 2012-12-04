function prep_modalsrc(osc, grid3d, eps_face_cell, mu_edge_cell, s_factor_cell, modalsrc)
chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(eps_face_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"eps_face_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(mu_edge_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"mu_edge_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count GK.count], [1 0]), ...
	'"s_factor_cell" should be %d-by-%d cell array whose each element is row vector with complex elements.', ...
	Axis.count, GK.count);
chkarg(istypesizeof(modalsrc, 'ModalSrc'), '"modalsrc" should be instance of ModalSrc.');
if ~isempty(modalsrc.Jh) && ~isempty(modalsrc.Jv)
	warning('FDS:srcAssign', 'Jh and Jv of "modalsrc" are already assigned; they will be reassigned.');
end

A = wgmode_matrix(osc.in_omega0(), ...
	eps_face_cell, mu_edge_cell, ...
	s_factor_cell, grid3d, ...
	modalsrc.normal_axis, modalsrc.intercept, ...
	PML.sc);
beta_guess = 2*pi*modalsrc.neff_guess/osc.in_L0();
[vec beta] = solve_mode(A, 'beta', beta_guess);
neff = beta*osc.in_L0()/2/pi;

vec = reshape(vec, 2, []).';
[h v] = cycle(modalsrc.normal_axis);
Hh = vec(:,1);
Hv = vec(:,2);
Nh = grid3d.N(h);
Nv = grid3d.N(v);

Hh = reshape(Hh, Nh, Nv);
Hv = reshape(Hv, Nh, Nv);

modalsrc.setJ(neff, Hv, -Hh, grid3d);
