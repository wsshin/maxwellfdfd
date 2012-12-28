function Dw = create_Dw(w, N, f1, fg)
% f1: factor multiplied to the first element in the w-direction
% fg: factor multiplied to the ghost element in the w-direction

% This function creates the forward derivative matrix.  For the backward
% derivative matrix, take the transpose or conjugate transpose appropriately.

chkarg(istypesizeof(w, 'Axis'), '"w" should be instance of Axis.');
chkarg(istypesizeof(N, 'int', [1 Axis.count]), '"N" should be length-%d row vector with integer elements.');
chkarg(istypesizeof(fg, 'complex'), '"nL" should be complex.');


shift = zeros(1, Axis.count);
shift(w) = 1;

% Get the displaced spatial markers.
ind_next = @(n, s) mod((1:n) + s - 1, n) + 1;
[i_next, j_next, k_next] = ndgrid(...
	ind_next(N(Axis.x), shift(Axis.x)), ...
	ind_next(N(Axis.y), shift(Axis.y)), ...
	ind_next(N(Axis.z), shift(Axis.z)));

% Translate spatial indices into matrix indices.
M = prod(N);
row_ind = 1:M;
col_ind_curr = 1:M;
col_ind_next = i_next + (j_next-1) * N(Axis.x) + (k_next-1) * N(Axis.x) * N(Axis.y);

a_curr = ones(N);
a_ind = {':', ':', ':'};
a_ind{w} = 1;
a_curr(a_ind{:}) = f1;

a_next = ones(N);
a_ind = {':', ':', ':'};
a_ind{w} = N(w);
a_next(a_ind{:}) = fg;

% Create the sparse matrix.
Dw = sparse([row_ind(:); row_ind(:)], [col_ind_curr(:); col_ind_next(:)], ...
            [-a_curr(:); a_next(:)], M, M);
