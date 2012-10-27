function [v, beta] = solve_mode(A, strategy, val, tol)

chkarg(isequal(strategy, 'beta') || isequal(strategy, 'modeindex'), '"strategy" should be ''beta'' or ''modeindex''.');

if nargin < 4  % no tol
	tol = 1e-6;
end

if isequal(strategy, 'beta')
	beta_guess = val;
	[v, lambda] = eigs(A, 1, -beta_guess^2);
	beta = sqrt(lambda)/1i;
else
	assert(isequal(strategy, 'modeindex'));
	mode_index = val;
	if isreal(A)
		[beta, v] = solve_mode_real(A, mode_index, tol, eps);
	else
		[beta_guess, ~] = solve_mode_real(A, mode_index, tol, 1e-2);
		[v, lambda] = eigs(A, 1, -real(beta_guess)^2);	
		beta = sqrt(-lambda);  % this guantees real(beta) >= 0
	end
%	[beta, v] = solve_mode_real2(A, mode_index);
end

function [beta, v] = solve_mode_real(A, mode_index, tol_extreme, tol_eigs)
% Assume that eigenvalues are all real, or at least quite close to real.  This
% is the case when the system has small loss.

chkarg(istypesizeof(mode_index, 'int') && mode_index >= 1, '"mode_index" should be positive integer.');

% Solve for largest-magnitude eigenvalue. 
%
% Find the largest magnitude of eigenvalues, and assume that this is the largest
% eigenvalue itself; if the frequency and loss of the system are small enough,
% the magnitude approximates the eigenvalue well.  We do this in order to solve
% for the most negative eigenvalues below, from which we can select the appropriate
% propagating mode.
lambda_max = max_eig_magnitude(A, tol_extreme);
lambda_max = lambda_max * 2;  % overestimate the maximum eigenvalue

% Shift the matrix and find the "mode_index" largest-magnitude eigenvalues and
% the corresponding eigenvectors.
B = A - lambda_max * speye(size(A,1));

lambda_min = -max_eig_magnitude(B, tol_extreme);  % note the (-) sign
if nargin < 4  % no tol_eigs
	opts.tol = eps;
else
	opts.tol = tol_eigs;
end
[V, D] = eigs(B, max(mode_index*2, 7), lambda_min, opts);  % calculate enough eigenvalues
% [V, D] = eigs(B, max(mode_index*2, 7), 'lm', opts);  % calculate enough eigenvalues
lambdas = diag(D);

betas = sqrt(-(lambdas + lambda_max));  % lambda = (i*beta)^2; betas have positive real parts.
[~, ind] = sort(real(betas), 'descend');

beta = betas(ind(mode_index));
v = V(:, ind(mode_index));

function [beta, v] = solve_mode_real2(A, mode_index)
% Assume that eigenvalues are all real, or at least quite close to real.  This
% is the case when the system has small loss.

chkarg(istypesizeof(mode_index, 'int') && mode_index >= 1, '"mode_index" should be positive integer.');

% Solve for largest-magnitude eigenvalue. 
%
% Find the largest magnitude of eigenvalues, and assume that this is the largest
% eigenvalue itself; if the frequency and loss of the system are small enough,
% the magnitude approximates the eigenvalue well.  We do this in order to solve
% for the most negative eigenvalues below, from which we can select the appropriate
% propagating mode.
%lambda_max = eigs(A, 1);
lambda_max = max_eig_magnitude(A, 1e-7);

% Shift the matrix and find the "mode_index" largest-magnitude eigenvalues and
% the corresponding eigenvectors.
B = A - lambda_max * speye(size(A,1));

opts.tol = 1e-14;
[V, D] = eigs(B, mode_index, 'lm', opts); 
lambdas = diag(D);
assert(all(lambdas <= 0));

beta = sqrt(lambdas(mode_index) + lambda_max)/1i;
v = V(:, mode_index);


function lambda_max = max_eig_magnitude(A, tol)
m = length(A);
if isreal(A)
	vmax = rand(m, 1);
else
	vmax = rand(2*m, 1);
	vmax = reshape(vmax, m, 2);
	vmax = vmax(:,1) + 1i * vmax(:,2);
end

lambda = norm(vmax);
relerr = inf;
while relerr >= tol
	vmax = vmax/lambda;
	vmax = A*vmax;
	lambda_old = lambda;
	lambda = norm(vmax);
	relerr = abs(lambda - lambda_old)/lambda_old;  % lambda_old > 0
end
lambda_max = lambda;
