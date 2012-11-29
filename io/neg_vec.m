function F_cell = neg_vec(F_cell)

chkarg(istypesizeof(F_cell, 'complexcell', [1 Axis.count]), ...
	'"F_cell" should be length-%d row cell array whose each element is array with complex elements.', Axis.count);

for w = Axis.elems
	F_cell{w} = -F_cell{w};
end
