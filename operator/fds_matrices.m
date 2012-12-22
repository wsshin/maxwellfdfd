%% fds_matrices
% Advanced feature which constructs the matrices and vectors that make up the simulation. 

%% Syntax
% [A1, A2, mu, epsilon, b] = fds_matrices(omega, s_prim, s_dual, mu, epsilon, J)

%% Description
% Produces the matrices and vectors needed to form the numerical model locally.
% This is an advanced feature that is useful if we want to locally check the
% error of the solution, or compute additional problem information such as
% derivative information.

%% Input parameters.
% Input parameters exactly match those used in the |fds.m| function.

%% Output parameters.
% * |A1|, |A2| -- matrices (sparse).
% * |m|, |e|, |b| -- vectors.

% Together, these allow one to form the problem on one's own machine.

%% Example: Local error check.
% Using the same parameters as from an |fds()| call, we can then compute the
% error in the simulation locally. Note that the locally computed error may
% differ from the server-computed value for implementation reasons.

%   % Run the simulation.
%   [E, H, err] = fds(omega, s_prim, s_dual, mu, epsilon, E, J, max_iters, err_thresh, 'text');
%
%   x = ([E{1}(:) ; E{2}(:) ; E{3}(:)]); % "Vectorize" E-field values.
%   y = ([H{1}(:) ; H{2}(:) ; H{3}(:)]); % "Vectorize" H-field values.
%
%   % Construct matrices
%   [A1, A2, m, e, b] = fds_matrices(omega, s_prim, s_dual, mu, epsilon, J); 
%
%   % Calculate errors in the wave equation for both E-fields and H-fields.
%   norm(A1 * ((1./m) .* (A2 * x)) - omega^2 * e .* x - b) / norm(b) % Error for E-field wave equation.
%   norm(A2 * ((1./e) .* (A1 * y)) - omega^2 * m .* y - A2 * (b ./ (-i*omega*e))) / norm(b) % Error for H-field wave equation.

function [A1, A2, m, e, b] = fds_matrices(omega, eps_face_cell, mu_edge_cell, s_factor_cell, J_cell, grid3d)

%% Get the shape.
N = grid3d.N;
M = prod(N);

s_prim = s_factor_cell(:, GK.prim);
s_dual = s_factor_cell(:, GK.dual);
d_prim = grid3d.dl(:, GK.prim);
d_dual = grid3d.dl(:, GK.dual);

%% Get the relevant derivative matrices.
[spx, spy, spz] = ndgrid(s_prim{Axis.x}.*d_prim{Axis.x}, s_prim{Axis.y}.*d_prim{Axis.y}, s_prim{Axis.z}.*d_prim{Axis.z});
[sdx, sdy, sdz] = ndgrid(s_dual{Axis.x}.*d_dual{Axis.x}, s_dual{Axis.y}.*d_dual{Axis.y}, s_dual{Axis.z}.*d_dual{Axis.z});

bc = grid3d.bc;

Df = cell(1, Axis.count);
for w = Axis.elems
	f1 = 1;
	
	if bc(w,Sign.p) == BC.p
		fg = exp(-1i * grid3d.kBloch(w) * grid3d.L(w));
	else  % bc(w,Sign.p) == BC.m
		fg = 0;
	end
	
	Df{w} = create_Dw(w, N, f1, fg);
end

Db = cell(1, Axis.count);
for w = Axis.elems
	if bc(w,Sign.n) == BC.e
		f1 = 2;
	else
		f1 = 1;
	end
	
	if bc(w,Sign.p) == BC.p
		fg = exp(-1i * grid3d.kBloch(w) * grid3d.L(w));
	else  % grid3d.bc(w,Sign.p) == BC.m
		fg = 0;
	end
	
	Db{w} = create_Dw(w, N, f1, fg);
	Db{w} = -Db{w}';  % conjugate transpose rather than transpose
end

Z = sparse(M, M);

my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));

% Forward differences (used to compute H from E).
Df{Axis.x} = my_diag(sdx.^-1) * Df{Axis.x};
Df{Axis.y} = my_diag(sdy.^-1) * Df{Axis.y};
Df{Axis.z} = my_diag(sdz.^-1) * Df{Axis.z};

% Backward differences (used to compute E from H).
Db{Axis.x} = my_diag(spx.^-1) * Db{Axis.x};
Db{Axis.y} = my_diag(spy.^-1) * Db{Axis.y};
Db{Axis.z} = my_diag(spz.^-1) * Db{Axis.z};

% Mask matrices
mask_m = cell(1, Axis.count);
for w = Axis.elems
	mask = ones(N);
	[u, v] = cycle(w);

	if bc(u,Sign.n) == BC.m
		ind = {':', ':', ':'};
		ind{u} = 1;
		mask(ind(:)) = 0;
	end
	if bc(v,Sign.n) == BC.m
		ind = {':', ':', ':'};
		ind{v} = 1;
		mask(ind(:)) = 0;
	end
	mask_m{w} = mask;
end
Mm = my_diag([mask_m{Axis.x}(:); mask_m{Axis.y}(:); mask_m{Axis.z}(:)]);

mask_e = cell(1, Axis.count);
for w = Axis.elems
	mask = ones(N);
	if bc(w,Sign.n) == BC.m
		ind = {':', ':', ':'};
		ind{w} = 1;
		mask(ind(:)) = 0;
	end
	mask_e{w} = mask;
end
Me = my_diag([mask_e{Axis.x}(:); mask_e{Axis.y}(:); mask_e{Axis.z}(:)]);

%% Form matrices
A1 = [Z, -Df{Axis.z}, Df{Axis.y}; ...
	Df{Axis.z}, Z, -Df{Axis.x}; ...
	-Df{Axis.y}, Df{Axis.x}, Z];
A1 = Me * A1 * Mm;

A2 = [Z, -Db{Axis.z}, Db{Axis.y}; ...
	Db{Axis.z}, Z, -Db{Axis.x}; ...
	-Db{Axis.y}, Db{Axis.x}, Z];
A2 = Mm * A2 * Me;

m = [mu_edge_cell{Axis.x}(:) ; mu_edge_cell{Axis.y}(:) ; mu_edge_cell{Axis.z}(:)];
e = [eps_face_cell{Axis.x}(:) ; eps_face_cell{Axis.y}(:) ; eps_face_cell{Axis.z}(:)];


b = -1i * omega * [J_cell{Axis.x}(:) ; J_cell{Axis.y}(:) ; J_cell{Axis.z}(:)];
