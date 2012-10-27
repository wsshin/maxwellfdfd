function prep_distsrc(osc, grid3d, eps_edge_cell, mu_face_cell, s_factor_cell, distsrc)
chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(eps_edge_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"eps_edge_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(mu_face_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"mu_face_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count GK.count], [1 0]), ...
	'"s_factor_cell" should be %d-by-%d cell array whose each element is row vector with complex elements.', ...
	Axis.count, GK.count);
chkarg(istypesizeof(distsrc, 'DistributedSrc'), '"distsrc" should be instance of DistributedSrc.');
if ~isempty(distsrc.Jh) && ~isempty(distsrc.Jv)
	warning('SrcAssign', 'Jh and Jv of "distsrc" are already assigned; they will be reassigned.');
end

A = wgmode_matrix(osc.in_omega0(), eps_edge_cell, mu_face_cell, s_factor_cell, grid3d, distsrc.normal_axis, distsrc.intercept, PML.SC);
beta_guess = 2*pi*distsrc.neff_guess/osc.in_L0();
[vec beta] = solve_mode(A, 'beta', beta_guess);
neff = beta*osc.in_L0()/2/pi;

vec = reshape(vec, 2, []).';
[h v] = cycle(distsrc.normal_axis);
Hh = vec(:,1);
Hv = vec(:,2);
Nh = grid3d.N(h);
Nv = grid3d.N(v);

Hh = reshape(Hh, Nh, Nv);
Hv = reshape(Hv, Nh, Nv);

distsrc.setJ(neff, Hv, -Hh, grid3d);
