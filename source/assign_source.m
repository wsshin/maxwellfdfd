function J_cell = assign_source(grid3d, source_array)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(source_array, 'Source', [1 0]), ...
	'"source_array" should be row vector of instances of Source.');

J_cell = cell(1, Axis.count);
for w = Axis.elems
	J_cell{w} = zeros(grid3d.N);
end

for src = source_array
	for w = Axis.elems
		[ind, Jw_patch] = src.generate(w, grid3d);
		if ~isempty(Jw_patch)
			J_cell{w}(ind{:}) = J_cell{w}(ind{:}) + Jw_patch;  % superpose sources
		end
	end
end
