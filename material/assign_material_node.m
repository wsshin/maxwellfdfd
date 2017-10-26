function [eps_node_cell, mu_node_cell] = assign_material_node(grid3d, object_array, eps_node_cell, mu_node_cell)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(object_array, 'EMObject', [1 0]), ...
	'"object_array" should be row vector of instances of EMObject.');

if nargin < 3  % no eps_node_cell
	eps_node_cell = {NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N)};
end
chkarg(istypesizeof(eps_node_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"eps_node_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.Ncell{:});

if nargin < 4  % no mu_node_cell
	mu_node_cell = {NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N)};
end
chkarg(istypesizeof(mu_node_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
	'"mu_node_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
	Axis.count, grid3d.Ncell{:});

ldual = cell(1, Axis.count);  % locations of cell centers
for w = Axis.elems
	ldual{w} = grid3d.l{w, GT.dual};  % grid3d.l rather than grid3d.lg
end

ind = cell(1, Axis.count);  % indices
for obj = object_array
	shape = obj.shape;
	material = obj.material;
	for w = Axis.elems
		bn = shape.bound(w, Sign.n);
		bp = shape.bound(w, Sign.p);
		in = find(ldual{w} >= bn, 1, 'first');
		ip = find(ldual{w} <= bp, 1, 'last');
		ind{w} = in:ip;
	end

	if istypesizeof(shape, 'Box')
		for w = Axis.elems
			eps_node_cell{w}(ind{:}) = material.eps(w,w);
			mu_node_cell{w}(ind{:}) = material.mu(w,w);
		end
	else  % shape is not a Box
		[X, Y, Z] = ndgrid(ldual{Axis.x}(ind{Axis.x}), ldual{Axis.y}(ind{Axis.y}), ldual{Axis.z}(ind{Axis.z}));
		is_in = shape.contains(X, Y, Z);
		ind_tf = false(grid3d.N);  % logical indices
		ind_tf(ind{:}) = is_in;
		
		for w = Axis.elems
			eps_node_cell{w}(ind_tf) = material.eps(w,w);
			mu_node_cell{w}(ind_tf) = material.mu(w,w);
		end
	end
end
