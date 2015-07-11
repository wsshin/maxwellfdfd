function Dw = create_Dw(w, N, f1, fg)
% f1: factor multiplied to the first element in the w-direction
% fg: factor multiplied to the ghost element in the w-direction

% This function creates the forward derivative matrix.  For the backward
% derivative matrix, take the transpose or conjugate transpose appropriately.

chkarg(istypesizeof(w, 'Axis') || istypesizeof(w, 'Dir'), '"w" should be instance of Axis or Dir.');
chkarg(istypesizeof(N, 'int', [1 w.count]), '"N" should be length-%d row vector with integer elements.', w.count);
chkarg(istypesizeof(fg, 'complex'), '"nL" should be complex.');

% Translate spatial indices (i,j,k) into matrix indices.
M = prod(N);
row_ind = 1:M;
col_ind_curr = 1:M;

col_ind_next = reshape(1:M, N);
shift = zeros(1, w.count);
shift(w) = -1;
col_ind_next = circshift(col_ind_next, shift);

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
