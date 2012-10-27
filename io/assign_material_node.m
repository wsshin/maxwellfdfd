function [eps_node_array, mu_node_array] = assign_material_node(grid3d, object_array)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(object_array, 'Object', [1 0]), ...
	'"object_array" should be row vector of instances of Object.');

eps_node_array = NaN(grid3d.N + 1);
mu_node_array = NaN(grid3d.N + 1);

lprim = cell(1, Axis.count);  % locations of the E-field grid planes
for w = Axis.elems
	lprim{w} = grid3d.lg{w, GK.prim};  % grid3d.lg rather than grid3d.l
end

ind = cell(1, Axis.count);  % indices
for obj = object_array
	shape = obj.shape;
	material = obj.material;
	for w = Axis.elems
		bn = shape.bound(w, Sign.n);
		bp = shape.bound(w, Sign.p);
		in = find(lprim{w} >= bn, 1, 'first');
		ip = find(lprim{w} <= bp, 1, 'last');
		ind{w} = in:ip;
	end

	if istypesizeof(shape, 'Box')
		eps_node_array(ind{Axis.x}, ind{Axis.y}, ind{Axis.z}) = material.eps;
		mu_node_array(ind{Axis.x}, ind{Axis.y}, ind{Axis.z}) = material.mu;
	else  % shape is not a Box
% 		for ix = ind{Axis.x}
% 			for iy = ind{Axis.y}
% 				for iz = ind{Axis.z}
% 					pt = [lprim{Axis.x}(ix), lprim{Axis.y}(iy), lprim{Axis.z}(iz)];
% 					if shape.contains(pt)
% 						eps_node(ix, iy, iz) = material.eps;
% 						mu_node(ix, iy, iz) = material.mu;
% 					end
% 				end
% 			end
% 		end
		[X, Y, Z] = ndgrid(lprim{Axis.x}(ind{Axis.x}), lprim{Axis.y}(ind{Axis.y}), lprim{Axis.z}(ind{Axis.z}));
		is_in = shape.contains([X(:), Y(:), Z(:)]);
		is_in = reshape(is_in, length(ind{Axis.x}), length(ind{Axis.y}), length(ind{Axis.z}));
		ind_tf = false(grid3d.N+1);  % logical indices
		ind_tf(ind{Axis.x}, ind{Axis.y}, ind{Axis.z}) = is_in;
		eps_node_array(ind_tf) = material.eps;
		mu_node_array(ind_tf) = material.mu;
	end
end
