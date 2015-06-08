% Stability test problem: $Projects/Trap and Release/maxwell/sidewall/symmetric thick DMD 1D/dispersion.m
% See June 4, 2015 note in $Projects/Trap and Release/Mathematica/RN - Trap and Release.nb 
% These should be included in the maxwellfdfd example directory in the
% future.

% Due to the random initialization of the guess eigenvector in eigs(), this
% function does not guarantee the consintent result.  For example, the
% above test example sometimes returns the wrong mode as the second-order
% mode.  However, tol_power = 1e-10 works pretty well.
function [v, beta] = solve_mode_real(A, modeorder, tol_power)
% Assume that eigenvalues are all real, or at least quite close to real.  This
% is the case when the system has small loss.

chkarg(ismatrix(A) && isreal(A), '"A" should be real matrix.')
chkarg(istypesizeof(modeorder, 'int') && modeorder > 0, '"modeorder" should be positive integer.');
chkarg(istypesizeof(tol_power, 'real') && tol_power > 0, '"tol_power" should be real and positive.');

% Solve for the largest-magnitude eigenvalue.
[~, lmax] = eig_largest_mag(A, tol_power);
assert(real(lmax) > 0);  % procedure below makes sense for lmax > 0
lmax = abs(lmax);  % make sure the shift quantity is real

% Eigenvalues are gamma = -beta^2, so they are negative.  We want to know
% the largest beta's, which correspond to the most negative gamma's.
% However, there is no way to directly calculate A's most negative
% eigenvalues by eigs(). Therefore, we shift A's eigenvalue spectrum by A's
% largest-magnitude eigenvalue. This makes all eigenvalues negative,
% because A's largest-magnitude eigenvalue is usually positive because the
% Nyquiest frequency is far larger than the magnitudes of the negative
% eigenvalues gamma (which are physical). Then, A's most negative
% eigenvalues become the largest-magnitude (and negative) eigenvalues,
% which can be calculated by eigs().
B = A - lmax * speye(length(A));  % all eigenvalues are negative

% The original code was to proceed to calculate [V, L] = eigs(B,
% modeorder). Because we are interested in A's most negative eigenvalues,
% and they correspond to B's largest-magnitude eigenvalues, the above code
% looked like a reasonable way to calculate the wanted eigenvalues.
% However, for some reason the above code often generated ARPACK errors.
% Forturately, I found a condition for calculating only one
% largest-magnitude eigenvalue of B stably.  The following procedure uses
% that eigenvalue of B.

% opts.v0 = rand(length(B),1);
% [~, lneg] = eigs(B, 1, 'lm', opts);  % lneg is most negative eigenvalue
[~, lneg] = eig_largest_mag(B, tol_power);  % lneg is most negative eigenvalue
lneg = -abs(lneg);  % make sure the shift quantity is real

% Below, lneg -> (1+sqrt(tol_power))*lneg makes sure that lshift shifts the
% eigenvalue spectrum of A to have all positive eigenvalues.  Too large
% factor (e.g., 2.00) generates the ARPACK error, because then the inverse
% iteration converges slowly.  (The inverse of the smallest and second
% smallest eigenvalues become very close.)
lshift = lmax + (1+sqrt(tol_power))*lneg;

B = A - lshift * speye(length(A));  % all eigenvalues are positive
[V, L] = eigs(B, modeorder+2, 'sm');  % two more eigenvalues just in case
ls = diag(L);

betas = sqrt(-(ls + lshift));
[~, ind] = sort(real(betas), 'descend');  % largest real part of beta first

beta = betas(modeorder);
v = V(:,ind(modeorder));
v = v / norm(v);


function [vmax, lmax] = eig_largest_mag(A, tol_power)
m = length(A);
if isreal(A)
	vmax = rand(m, 1);
else
	vmax = rand(2*m, 1);
	vmax = reshape(vmax, m, 2);
	vmax = vmax(:,1) + 1i * vmax(:,2);
end
vmax = vmax / norm(vmax);  % vmax is normalized

% Calculate the maximum eigenvalue by power iteration.
lmax_prev = Inf;
vmax_next = A * vmax;
lmax = vmax' * vmax_next;  % Rayleigh quotient
vmax_next = vmax_next / norm(vmax_next);  % vmax_next is normalized
while abs((lmax - lmax_prev) / lmax) > tol_power
	vmax = vmax_next;
	lmax_prev = lmax;
	vmax_next = A * vmax;
	lmax = vmax' * vmax_next;
	vmax_next = vmax_next / norm(vmax_next);  % vmax_next is normalized
end
