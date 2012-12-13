function Prod_cell = mult_vec(F_cell, G_cell)

chkarg(istypesizeof(F_cell, 'complexcell', [1 Axis.count]) || istypesizeof(F_cell, 'complexcell', [Axis.count 1]), ...
	'"F_cell" should be length-%d row cell array whose each element is array with complex elements.', Axis.count);
chkarg(istypesizeof(G_cell, 'complexcell', [1 Axis.count]) || istypesizeof(G_cell, 'complexcell', [Axis.count 1]), ...
	'"F_cell" should be length-%d row cell array whose each element is array with complex elements.', Axis.count);

Prod_cell = cell(1, Axis.count);
for w = Axis.elems
	chkarg(all(size(F_cell{w})==size(G_cell{w})), '%d-component of "F_cell" and "G_cell" should have the same size.');
	Prod_cell{w} = F_cell{w} .* G_cell{w};
end
