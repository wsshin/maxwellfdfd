function F_cell = flip_vec(F_cell)

chkarg(istypesizeof(F_cell, 'complexcell', [1 Axis.count]) || istypesizeof(F_cell, 'complexcell', [Axis.count 1]), ...
	'"F_cell" should be length-%d row or column cell array whose each element is array with complex elements.', Axis.count);

for w = Axis.elems
	F_cell{w} = flip_array(F_cell{w});
end
