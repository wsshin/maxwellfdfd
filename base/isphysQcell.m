function truth = isphysQcell(cell_array)
truth = iscell(cell_array) && ismatrix(cell_array) && size(cell_array, 2) == 2;

n = size(cell_array, 1);
for i = 1:n
	truth = truth && istypesizeof(cell_array{i,1}, 'PhysQ') && istypesizeof(cell_array{i,2}, 'int');
end
