function [E, H, A, b, HfromE] = solve_eq_direct(omega, eps_face_cell, mu_edge_cell, s_factor_cell, J_cell, grid3d, nosolve)

if nargin < 7  % no nosolve
	nosolve = false;
end

[A1, A2, mu, eps, b] = fds_matrices(omega, eps_face_cell, mu_edge_cell, s_factor_cell, J_cell, grid3d);

% Mask elements corresponding to PEC.
ind_pec = isinf(abs(eps));
eps(ind_pec) = 1;
pec_mask = ones(size(ind_pec));
pec_mask(ind_pec) = 0;
PM = spdiags(pec_mask, 0, length(pec_mask), length(pec_mask));

INV_MU = spdiags(1./mu, 0, length(mu), length(mu));  % when mu has Inf, "MU \ Mat" complains about singularity
EPS = spdiags(eps, 0, length(eps), length(eps));

HfromE = (INV_MU * A2);
A = PM * A1 * HfromE * PM - omega^2 * EPS;
HfromE = (1i/omega) * HfromE;

% Reorder the indices of the elements of matrices and vectors to reduce the bandwidth of A.
N = grid3d.N;
r = reordering_indices(Axis.count, N);

A = A(r,r);
b = b(r);
HfromE = HfromE(r,r);

% figure; spy(A); xlabel(''); set(gca, 'xtick', []); set(gca, 'ytick', []);
% [L1, U1] = lu(A);
% figure; spy(U1); xlabel(''); set(gca, 'xtick', []); set(gca, 'ytick', []);
% [L2, U2, p, q, R] = lu(A, 'vector');
% figure; spy(U2); xlabel(''); set(gca, 'xtick', []); set(gca, 'ytick', []);


if ~nosolve
	% Below, [L, U, P, Q, R] = lu(A), rather than [L, U] = lu(A), seems to be
	% similar to this.  The permutation matrices P and Q can be generated as
	% vectors by [L, U, p, q, R] = lu(A, 'vector').  Then R(:,p)\A(:,q) = L*U.
	
%	spparms('spumoni',1);  % make mldivide(A, b) (= A\b) verbose
	e = A\b;  
	h = HfromE * e;

	e = reshape(e, Axis.count, prod(N));
	Ex = e(int(Axis.x), :); Ex = reshape(Ex, N);
	Ey = e(int(Axis.y), :); Ey = reshape(Ey, N); 
	Ez = e(int(Axis.z), :); Ez = reshape(Ez, N);
	E = {Ex, Ey, Ez};

	% Test symmetry with respect to the plane bisecting the x-axis.
	if false
		E = test_sym(Axis.x, A, b, E);
	end

	h = reshape(h, Axis.count, prod(N));
	Hx = h(int(Axis.x), :); Hx = reshape(Hx, N);
	Hy = h(int(Axis.y), :); Hy = reshape(Hy, N); 
	Hz = h(int(Axis.z), :); Hz = reshape(Hz, N);
	H = {Hx, Hy, Hz};
else
	E = {};
	H = {};
end

function E = test_sym(w, A, b, E)
Nw = size(E{Axis.x}, int(w));

Ex = E{Axis.x}; Ey = E{Axis.y}; Ez = E{Axis.z};
e = [Ex(:), Ey(:), Ez(:)];
e = e.';
e = e(:);
fprintf('norm(b - A*e)/norm(b) = %e\n', norm(b - A*e)/norm(b));

if mod(Nw, 2) == 0
	ind2 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind2{v}{w} = [1:Nw/2+1, Nw/2:-1:2];
		else
			ind2{v}{w} = [1:Nw/2, Nw/2:-1:1];
		end
	end
	
	ind3 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind3{v}{w} = [1, Nw:-1:Nw/2+1, Nw/2+2:Nw];
		else
			ind3{v}{w} = [Nw:-1:Nw/2+1, Nw/2+1:Nw];
		end
	end
else
	Nw_half = (Nw+1)/2;
	
	ind2 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind2{v}{w} = [1:Nw_half, Nw_half:-1:2];
		else
			ind2{v}{w} = [1:Nw_half, Nw_half-1:-1:1];
		end
	end
	
	ind3 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind3{v}{w} = [1, Nw:-1:Nw_half+1, Nw_half+1:Nw];
		else
			ind3{v}{w} = [Nw:-1:Nw_half, Nw_half+1:Nw];
		end
	end
end
Ex2 = Ex(ind2{Axis.x}{:});
Ey2 = Ey(ind2{Axis.y}{:});
Ez2 = Ez(ind2{Axis.z}{:});
E2 = {Ex2, Ey2, Ez2};

Ex3 = Ex(ind3{Axis.x}{:});
Ey3 = Ey(ind3{Axis.y}{:});
Ez3 = Ez(ind3{Axis.z}{:});
E3 = {Ex3, Ey3, Ez3};

Ex4 = (Ex2 + Ex3)/2;
Ey4 = (Ey2 + Ey3)/2;
Ez4 = (Ez2 + Ez3)/2;
E4 = {Ex4, Ey4, Ez4};

e2 = [Ex2(:), Ey2(:), Ez2(:)];
e2 = e2.';
e2 = e2(:);

e3 = [Ex3(:), Ey3(:), Ez3(:)];
e3 = e3.';
e3 = e3(:);

e4 = [Ex4(:), Ey4(:), Ez4(:)];
e4 = e4.';
e4 = e4(:);

fprintf('norm(b - A*e2)/norm(b) = %e\n', norm(b - A*e2)/norm(b));
fprintf('norm(b - A*e3)/norm(b) = %e\n', norm(b - A*e3)/norm(b));
fprintf('norm(b - A*e4)/norm(b) = %e\n', norm(b - A*e4)/norm(b));
fprintf('norm(b) = %e\n', norm(b));
fprintf('norm(e) = %e\n', norm(e));
fprintf('norm(A, 1) = %e\n', norm(A, 1));

E = E3;
