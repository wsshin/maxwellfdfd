function prep_modalsrc(ge, pml, osc, grid3d, eps_cell, mu_cell, s_factor_cell, modalsrc)
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');
chkarg(istypesizeof(pml, 'PML'), '"pml" should be instance of PML.');
chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(eps_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"eps_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(mu_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"mu_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
	'"s_factor_cell" should be %d-by-%d cell array whose each element is row vector with complex elements.', ...
	Axis.count, GT.count);
chkarg(istypesizeof(modalsrc, 'ModalSrc'), '"modalsrc" should be instance of ModalSrc.');
if ~isempty(modalsrc.JMh) && ~isempty(modalsrc.JMv)
	warning('FDS:srcAssign', 'Jh and Jv of "modalsrc" are already assigned; they will be reassigned.');
end

[A, En_Ht, gEh_HvEn, gEv_HhEn, Hn_Et] = wgmode_matrix(ge, pml, osc.in_omega0(), ...
	eps_cell, mu_cell, ...
	s_factor_cell, grid3d, ...
	modalsrc.normal_axis, modalsrc.intercept);
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

modalsrc.setEH(neff, osc, E_cell, H_cell, ge, grid3d);
