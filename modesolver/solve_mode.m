% This function is no longer used; see prep_modalsrc(). 
function [v, beta] = solve_mode(A, strategy, sval, v_guess)

chkarg(isequal(strategy, 'beta') || isequal(strategy, 'modeorder'), '"strategy" should be ''beta'' or ''modeorder''.');
if isequal(strategy, 'beta')
	chkarg(istypesizeof(sval, 'complex'), '"sval" should be "beta" (complex).');
else  % strategy == 'modeorder'
	chkarg(istypesizeof(sval, 'int') && sval > 0, '"sval" should be "modeorder" (positive integer).');
end

if isequal(strategy, 'beta')
	beta_guess = sval;
	
	% eigs() generates an error if the shifte matrix is sigular, i.e., if the
	% shift value is actually an eigenvalue of A.  To prevent such an error,
	% perturb the shift value by adding eps.
	[v, l] = eigs(A, 1, -beta_guess^2 + eps);  % gamma^2 = (i beta)^2
	beta = sqrt(-l);  % this guantees real(beta) >= 0
else  % strategy == modeorder
	modeorder = sval;
	assert(isreal(A));
	[beta, v] = solve_mode_real(A, modeorder, 1e-10);
	
% 	if isreal(A)
% 		[beta, v] = solve_mode_real(A, modeorder, 1e-10);
% 	else
% 		[beta_guess, v] = solve_mode_real(real(A), modeorder, 1e-4);
% 		opts.v0 = v_guess;
% 		[v, l] = eigs(A, 1, -real(beta_guess)^2, opts);	
% 
% 		% Perform the Rayleight quotient iteration on A with x as an initial guess for
% 		% the deiserd eigenvector.  (This process can be justified only for sufficiently
% 		% small loss.)
% 		l_prev = -real(beta_guess)^2;
% 		l = v' * A * v;
% 		i = 0;
% 		n = length(A);
% 		tol = 1e-13;
% 		while abs((l - l_prev) / l) > tol && i < 20
% 			fprintf('%s\n', num2str(abs((l - l_prev) / l)));
% 			i = i + 1;
% 			l_prev = l;
% 			v = (A - l * speye(n)) \ v;
% 			v = v / norm(v);  % v is normalized
% 			l = v' * A * v;
% 		end
		beta = sqrt(-l);  % this guantees real(beta) >= 0
	end
end

function [beta, v] = solve_mode_real(A, modeorder, tol_power)
% Assume that eigenvalues are all real, or at least quite close to real.  This
% is the case when the system has small loss.

chkarg(istypesizeof(modeorder, 'int') && modeorder > 0, '"modeorder" should be positive integer.');

% Solve for the largest-magnitude eigenvalue.
lmax = eig_largest_mag(A, tol_power);

% Shift the matrix and find the "mode_index" largest-magnitude eigenvalues and
% the corresponding eigenvectors.
B = A - lmax * speye(length(A));

if modeorder == 1
	[l, v] = eig_largest_mag(B, tol_power);
else  % modeorder >= 2
	[V, L] = eigs(B, modeorder);
	ls = diag(L);
	[~, ind] = sort(ls);  % most negative eigenvalues first.

	l = ls(modeorder);
	v = V(:,ind(modeorder));
end

beta = sqrt(-(l + lmax));
v = v / norm(v);


function [lmax, vmax] = eig_largest_mag(A, tol_power)
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
