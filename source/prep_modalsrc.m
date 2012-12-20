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

[A, En_Ht, gEh_HvEn, gEv_HhEn, Hn_Et] = wgmode_matrix(osc.in_omega0(), ...
	eps_face_cell, mu_edge_cell, ...
	s_factor_cell, grid3d, ...
	modalsrc.normal_axis, modalsrc.intercept, ...
	PML.sc);
beta_guess = 2*pi*modalsrc.neff_guess/osc.in_L0();
[Ht, beta] = solve_mode(A, 'beta', beta_guess);
neff = beta*osc.in_L0()/2/pi;

Ht = reshape(Ht, Dir.count, []).';
[h, v, n] = cycle(modalsrc.normal_axis);
Hh = Ht(:, Dir.h);
Hv = Ht(:, Dir.v);

En = En_Ht * [Hh; Hv];
Eh = gEh_HvEn * [Hv; En] ./ (1i*beta);  % gamma = 1i * beta
Ev = gEv_HhEn * [Hh; En] ./ (1i*beta);  % gamma = 1i * beta
Hn = Hn_Et * [Eh; Ev];

Nh = grid3d.N(h);
Nv = grid3d.N(v);

Eh = reshape(Eh, Nh, Nv);
Ev = reshape(Ev, Nh, Nv);
En = reshape(En, Nh, Nv);

Hh = reshape(Hh, Nh, Nv);
Hv = reshape(Hv, Nh, Nv);
Hn = reshape(Hn, Nh, Nv);

E_cell = cell(1, Axis.count);
E_cell{h} = Eh;
E_cell{v} = Ev;
E_cell{n} = En;

H_cell = cell(1, Axis.count);
H_cell{h} = Hh;
H_cell{v} = Hv;
H_cell{n} = Hn;

modalsrc.setEH(neff, osc, E_cell, H_cell, grid3d);
