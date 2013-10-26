function JM_cell = assign_source(grid3d, src_array)


chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(src_array, 'Source', [1 0]), ...
	'"src_array" should be row vector of instances of Source.');

JM_cell = cell(1, Axis.count);
for w = Axis.elems
	JM_cell{w} = zeros(grid3d.N);
end

for src = src_array
	for w = Axis.elems
		[ind, JMw_patch] = src.generate(w, grid3d);
		if ~isempty(JMw_patch)
			JM_cell{w}(ind{:}) = JM_cell{w}(ind{:}) + JMw_patch;  % superpose sources
		end
	end
end
