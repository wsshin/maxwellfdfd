function [E, H, A, x0, b, HfromE] = solve_eq_direct(osc, grid3d, s_factor, mu_face, eps_edge, J, E0, nosolve)

if nargin < 8  % no nosolve
	nosolve = false;
end

[A1, A2, mu, eps, b] = fds_matrices(osc.in_omega0(), ...
						grid3d.dl(:,GK.prim), grid3d.dl(:,GK.dual), ...
						s_factor(:,GK.prim), s_factor(:,GK.dual), ...
						mu_face, eps_edge, ...
						J);

% Mask elements corresponding to PEC.
ind_pec = isinf(abs(eps));
eps(ind_pec) = 1;
pec_mask = ones(size(ind_pec));
pec_mask(ind_pec) = 0;
PM = spdiags(pec_mask, 0, length(pec_mask), length(pec_mask));

MU = spdiags(mu, 0, length(mu), length(mu));
EPS = spdiags(eps, 0, length(eps), length(eps));
omega = osc.in_omega0();

HfromE = (MU \ A2);
A = PM * A1 * HfromE * PM - omega^2 * EPS;
HfromE = (1i/omega) * HfromE;
x0 = [E0{Axis.x}(:); E0{Axis.y}(:); E0{Axis.z}(:)];

if ~nosolve
	e = A\b;
	h = HfromE * e;

	e = reshape(e, length(e)/3, []);
	Ex = e(:, 1); Ex = reshape(Ex, grid3d.N);
	Ey = e(:, 2); Ey = reshape(Ey, grid3d.N); 
	Ez = e(:, 3); Ez = reshape(Ez, grid3d.N);
	E = {Ex, Ey, Ez};

	h = reshape(h, length(h)/3, []);
	Hx = h(:, 1); Hx = reshape(Hx, grid3d.N);
	Hy = h(:, 2); Hy = reshape(Hy, grid3d.N); 
	Hz = h(:, 3); Hz = reshape(Hz, grid3d.N);
	H = {Hx, Hy, Hz};
else
	E = {};
	H = {};
end
