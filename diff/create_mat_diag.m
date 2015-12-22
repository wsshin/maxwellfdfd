function mat_diag_cell = create_mat_diag(mat_array, grid3d)

chkarg(istypesizeof(mat_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
	'"mat_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});

mat_diag_cell = cell(1, Axis.count);

Ntot = prod(grid3d.N);  % total number of grid cells

diag_vec = NaN(Ntot, 1);
for w = Axis.elems
	k = int(w) - 1;  % kth off-diagonal
	for v = Axis.elems
		r = int(v);  % row index
		c = mod(r + k - 1, Axis.count) + 1;  % column index
		diag_array = mat_array(r, c, :, :, :);
		diag_vec(1+(r-1)*Ntot : r*Ntot) = diag_array(:);
	end
	
	mat_diag_cell{w} = create_spdiag(diag_vec);
end
