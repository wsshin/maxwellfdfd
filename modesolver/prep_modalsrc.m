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
	warning('Maxwell:srcAssign', 'Jh and Jv of "modalsrc" are already assigned; they will be reassigned.');
end

% % Assuming that loss is small, construct the waveguide mode matrix ignoring the
% % lossy components.
% Ar = wgmode_matrix(ge, pml, real(osc.in_omega0()), ...
% 	real_cell(eps_cell), real_cell(mu_cell), ...
% 	real_cell(s_factor_cell), grid3d, ...
% 	modalsrc.normal_axis, modalsrc.intercept);
% 
% % Ar is not symmetric on a nonuniform grid.  However, we can symmetrize it using
% % a similarity transformation with a real diagonal matrix.  The symmetrized
% % matrix is, therefore, still real, and hence has all-real eigenvalues.  Because
% % similar matrices have the same eigenvalues, Ar should also have all-real
% % eigenvalues, and therefore shifting the eigenvalue spectrum below makes sense.
% 
% % Calculate the maximum eigenvalue by power iteration.
% n = size(Ar, 1);
% lmax_prev = Inf;
% xmax = randn(n, 1);
% xmax = xmax / norm(xmax);  % vmax is normalized
% xmax_prev = xmax;
% xmax = Ar * xmax;
% lmax = xmax_prev' * xmax;  % Rayleigh quotient
% xmax = xmax / norm(xmax);  % vmax is normalized
% i = 0;
% tol = 1e-4;
% while abs((lmax - lmax_prev) / lmax) > tol && i < 100
% 	i = i + 1;
% 	lmax_prev = lmax;
% 	xmax_prev = xmax;
% 	xmax = Ar * xmax;
% 	lmax = xmax_prev' * xmax;
% 	xmax = xmax / norm(xmax);  % xmax is normalized
% end
% 
% if abs((lmax - lmax_prev) / lmax) > tol
% 	warning('Maxwell:modeSolver', 'power iteration did not converge: abs((lmax - lmax_prev) / lmax) = %s', ...
% 		num2str(abs((lmax - lmax_prev) / lmax)));
% end
% 
% % shift = abs(lmax);
% shift = lmax;
% modeorder = modalsrc.sval;
% [X, L] = eigs(Ar - shift * speye(n), modeorder);  % calculate a few more modes to be sure
% ls = diag(L);
% [~, ind] = sort(ls);  % most negative eigenvalues first.
% l = ls(modeorder) + shift;
% x = X(:,ind(modeorder));
% x = x / norm(x);  % x is normalized
% 
% 
% % Construct the waveguide mode matrix including lossy components.
% [A, En_Htr, gEh_HvEn, gEv_HhEn, Hn_Etr] = wgmode_matrix(ge, pml, osc.in_omega0(), ...
% 	eps_cell, mu_cell, ...
% 	s_factor_cell, grid3d, ...
% 	modalsrc.normal_axis, modalsrc.intercept);
% 
% % Perform the Rayleight quotient iteration on A with x as an initial guess for
% % the deiserd eigenvector.  (This process can be justified only for sufficiently
% % small loss.)
% l_prev = l;
% l = x' * A * x;
% i = 0;
% tol = 1e-13;
% while abs((l - l_prev) / l) > tol && i < 20
% 	fprintf('%s\n', num2str(abs((l - l_prev) / l)));
% 	i = i + 1;
% 	l_prev = l;
% 	x = (A - l * speye(n)) \ x;
% 	x = x / norm(x);  % x is normalized
% 	l = x' * A * x;
% end
% 
% if abs((l - l_prev) / l) > tol
% warning('Maxwell:modeSolver', 'Rayleigh quotient iteration did not converge: abs((l - l_prev) / l) = %s', ...
% 	num2str(abs((l - l_prev) / l)));
% end
% 
% Htr = x;
% beta = sqrt(-l);  % this guarantees real(beta) >= 0; note that l = gamma^2 = -beta^2
% neff = beta*osc.in_L0()/2/pi;

[A, En_Htr, gEh_HvEn, gEv_HhEn, Hn_Etr] = wgmode_matrix(ge, pml, osc.in_omega0(), ...
	eps_cell, mu_cell, ...
	s_factor_cell, grid3d, ...
	modalsrc.normal_axis, modalsrc.intercept);

if isequal(modalsrc.opts.clue, 'order')
	Ar = wgmode_matrix(ge, pml, real(osc.in_omega0()), ...
		real_cell(eps_cell), real_cell(mu_cell), ...
		real_cell(s_factor_cell), grid3d, ...
		modalsrc.normal_axis, modalsrc.intercept);

	[v, beta_guess] = solve_mode_real(Ar, modalsrc.opts.order, 1e-6);
	opts.v0 = v;
else
	assert(isequal(modalsrc.opts.clue, 'guess'))
	beta_guess = 2*pi*modalsrc.opts.neff/osc.in_L0();
	if isfield(modalsrc.opts, 'H2d')
		grid2d = Grid2d(grid3d, modalsrc.normal_axis);
		Hi = cell(1, Dir.count);  % interpreted H field
		for d = Dir.elems
			w = grid2d.axis(d);
			Hw2d = modalsrc.opts.H2d{w};
			
			[Hw, l] = Hw2d.data_expanded();
			[Xh, Yv] = ndgrid(l{:});
			
			li = grid2d.l(Dir.elems + Dir.count*subsindex(Hw2d.gt_array));
			[Xhi, Yvi] = ndgrid(li{:});
			Hi{d} = interpn(Xh, Yv, Hw, Xhi, Yvi);
		end
		v = [Hi{Dir.h}(:), Hi{Dir.v}(:)].';
		v = v(:);
		opts.v0 = v;
	end
end

opts.tol = eps;  % MATLAB default value; just to make sure opts exists

use_eigs = false;
if ~isfield(modalsrc.opts, 'H2d')  || use_eigs  % no guess on field distribution
	% eigs() generates an error if the shift matrix is sigular, i.e., if the shift
	% value is actually an eigenvalue of A.  To prevent such an error, perturb the
	% shift value by adding eps.
	[v, l] = eigs(A, 1, -beta_guess^2 + eps, opts);  % gamma^2 = (i beta)^2
else  % there is guess on field distribution
	% Perform the Rayleight quotient iteration on A with x as an initial guess for
	% the desired eigenvector.  (This process can be justified only for sufficiently
	% small loss, because beta_guess is for the lossless system.)
	l_prev = -beta_guess^2;
	l = v' * A * v;
	i = 0;
	tol = 1e-13;
	n = length(A);
	while abs((l - l_prev) / l) > tol && i < 20
% 		fprintf('\t%s\n', num2str(abs((l - l_prev) / l)));
		i = i + 1;
		l_prev = l;
		v = (A - l * speye(n)) \ v;
		v = v / norm(v);  % x is normalized
		l = v' * A * v;
	end
end

beta = sqrt(-l);
Htr = v;
neff = beta*osc.in_L0()/2/pi;

Htr = reshape(Htr, Dir.count, []).';
[h, v, n] = cycle(modalsrc.normal_axis);
Hh = Htr(:, Dir.h);
Hv = Htr(:, Dir.v);

%% En from Htr
En = En_Htr * [Hh; Hv];

%% Etr from Htr and En
Eh = gEh_HvEn * [Hv; En] ./ (1i*beta);  % gamma = 1i * beta
Ev = gEv_HhEn * [Hh; En] ./ (1i*beta);  % gamma = 1i * beta

%% Hn from Etr
Hn = Hn_Etr * [Eh; Ev];

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
