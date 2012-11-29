function r = reordering_indices(dof, N)
% Generate indices to reorder the elements of matrices and vectors to reduce the
% bandwidth of the Maxwell operator matrix.

chkarg(istypesizeof(dof, 'int') && dof > 0, '"dof" should be positive integer.');
chkarg(istypesizeof(N, 'int', [1 0]), '"N" should be row vector with integer elements.');

r = 1:dof*prod(N);
r = reshape(r, prod(N), dof);
r = r.';
r = r(:);
