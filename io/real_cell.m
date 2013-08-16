function real_cell = real_cell(complex_cell)

chkarg(istypeof(complex_cell, 'complexcell'), ...
	'"complex_cell" should be cell array with complex elements.');

real_cell = cell(size(complex_cell));
n = length(real_cell(:));
for i = 1:n
	real_cell{i} = real(complex_cell{i});
end
