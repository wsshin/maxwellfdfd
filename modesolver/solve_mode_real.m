function [v, beta] = solve_mode_real(A, modeorder, tol_power)
% Assume that eigenvalues are all real, or at least quite close to real.  This
% is the case when the system has small loss.

chkarg(ismatrix(A) && isreal(A), '"A" should be real matrix.')
chkarg(istypesizeof(modeorder, 'int') && modeorder > 0, '"modeorder" should be positive integer.');
chkarg(istypesizeof(tol_power, 'real') && tol_power > 0, '"tol_power" should be real and positive.');

% Solve for the largest-magnitude eigenvalue.
[~, lmax] = eig_largest_mag(A, tol_power);

% Shift the matrix and find the "mode_index" largest-magnitude eigenvalues and
% the corresponding eigenvectors.
B = A - lmax * speye(length(A));

[V, L] = eigs(B, modeorder);  % calculating more modes with (modeorder + 2) does not work for single-mode waveguide
ls = diag(L);
[~, ind] = sort(ls);  % most negative eigenvalues first.

l = ls(modeorder);
v = V(:,ind(modeorder));

beta = sqrt(-(l + lmax));
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
