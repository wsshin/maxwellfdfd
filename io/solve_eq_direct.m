function [E, H, A, x0, b, HfromE] = solve_eq_direct(omega, d_prim, d_dual, s_prim, s_dual, mu_face, eps_edge, J, E0, nosolve)

if nargin < 10  % no nosolve
	nosolve = false;
end

N = [length(d_prim{Axis.x}), length(d_prim{Axis.y}), length(d_prim{Axis.z})];

[A1, A2, mu, eps, b] = fds_matrices(omega, ...
						d_prim, d_dual, ...
						s_prim, s_dual, ...
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

HfromE = (MU \ A2);
A = PM * A1 * HfromE * PM - omega^2 * EPS;
HfromE = (1i/omega) * HfromE;
x0 = [E0{Axis.x}(:); E0{Axis.y}(:); E0{Axis.z}(:)];

if ~nosolve
	e = A\b;
	h = HfromE * e;

	e = reshape(e, length(e)/3, []);
	Ex = e(:, 1); Ex = reshape(Ex, N);
	Ey = e(:, 2); Ey = reshape(Ey, N); 
	Ez = e(:, 3); Ez = reshape(Ez, N);
	E = {Ex, Ey, Ez};

	h = reshape(h, length(h)/3, []);
	Hx = h(:, 1); Hx = reshape(Hx, N);
	Hy = h(:, 2); Hy = reshape(Hy, N); 
	Hz = h(:, 3); Hz = reshape(Hz, N);
	H = {Hx, Hy, Hz};
else
	E = {};
	H = {};
end
