function Prod_cell = mult_vec(F_cell, G_cell)

chkarg(istypesizeof(F_cell, 'complexcell', [1 0]) || istypesizeof(F_cell, 'complexcell', [0 1]), ...
	'"F_cell" should be row or column cell array whose each element is array with complex elements.');
chkarg(istypesizeof(G_cell, 'complexcell', [1 0]) || istypesizeof(G_cell, 'complexcell', [0 1]), ...
	'"F_cell" should be row or column cell array whose each element is array with complex elements.', Axis.count);
chkarg(length(F_cell)==length(G_cell), '"F_cell" and "G_cell" should be the same length.');

n = length(F_cell);

Prod_cell = cell(size(F_cell));
for w = 1:n
	chkarg(all(size(F_cell{w})==size(G_cell{w})), '%d-component of "F_cell" and "G_cell" should have the same size.');
	Prod_cell{w} = F_cell{w} .* G_cell{w};
end
