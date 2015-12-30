function [eps_array, mu_array] = mean_subpixel(grid3d, ge, domain, eps_ind_cell, mu_ind_cell, eps_shape_ind_cell, mu_shape_ind_cell, ind2eps_array, ind2mu_array, ind2shape_array, eps_array, mu_array)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');
chkarg(istypesizeof(domain, 'Domain'), '"domain" shuold be instance of Domain.');
chkarg(istypesizeof(eps_ind_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"eps_ind_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(mu_ind_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"mu_ind_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(eps_shape_ind_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"eps_shape_ind_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(mu_shape_ind_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"mu_shape_ind_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer as elements.', ...
	Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(ind2eps_array, 'complex', [Axis.count Axis.count 0]), ...
	'"ind2eps_array" should be %d-by-%d-by-n array with complex elements.', Axis.count, Axis.count);
chkarg(istypesizeof(ind2mu_array, 'complex', [Axis.count Axis.count 0]), ...
	'"ind2mu_array" should be %d-by-%d-by-n array with complex elements.', Axis.count, Axis.count);
chkarg(istypesizeof(ind2shape_array, 'Shape', [1 0]), ...
	'"ind2shape_array" should be row vector with Shape as elements.');
chkarg(istypesizeof(eps_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
	'"eps_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(mu_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
	'"mu_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});

% Permute eps_array and mu_array to use 2D logical indices.  Logical indices of
% more than 1D can be used only as the last indices.
eps_array = permute(eps_array, [3 4 5 1 2]);
mu_array = permute(mu_array, [3 4 5 1 2]);

[I, J, K] = ndgrid(1:grid3d.N(Axis.x), 1:grid3d.N(Axis.y), 1:grid3d.N(Axis.z));
union_shape_map = containers.Map();
for ft = FT.elems  % eps or mu
	if ft == FT.e  % eps
		gt = alter(ge);  % for smoothing eps, voxels to examine have H-fields at vertices
		mat_ind_cell = eps_ind_cell;
		shape_ind_cell = eps_shape_ind_cell;
		ind2mat_array = ind2eps_array;
	else  % mu
		gt = ge;  % for smoothing mu, voxels to examine have E-fields at vertices
		mat_ind_cell = mu_ind_cell;
		shape_ind_cell = mu_shape_ind_cell;
		ind2mat_array = ind2mu_array;
	end
	
	% Calculate the diagonal entries of eps and mu tensors.
	for w = [1+Axis.count int(Axis.elems)]
		gt_array = gt(ones(1, Axis.count));
		if w ~= 1+Axis.count  % w == x, y, z
			gt_array(w) = alter(gt);
		end
		
		mat_ind_array = mat_ind_cell{w};
		shape_ind_array = shape_ind_cell{w};
		
		mat_ind_array = attach_ghost_mat_ind(mat_ind_array, gt_array, grid3d.bc);
		[shape_ind_array, ind2shape_array] = attach_ghost_shape_ind(domain, shape_ind_array, ind2shape_array, gt_array, grid3d.bc);

		is_hetero_mat_array = find_hetero_voxels(grid3d, mat_ind_array);  % imagesc(is_hetero_mat_array) shows locations to smooth

		l_voxel = grid3d.lg(Axis.elems + Axis.count*subsindex(gt_array));  % locations of vertices of voxels
		l_voxel_center = grid3d.l(Axis.elems + Axis.count*subsindex(alter(gt_array)));  % locations of quasi-centers of voxels; note this is l, not lg (consider primary voxel)
		
		loc_ind_hetero_array = [I(is_hetero_mat_array), J(is_hetero_mat_array), K(is_hetero_mat_array)];  % n_hetero-by-3
		n_hetero = size(loc_ind_hetero_array, 1);  % number of hetero voxels
		
		for m = 1:n_hetero
			loc_ind = num2cell(loc_ind_hetero_array(m, :));
			[i, j, k] = deal(loc_ind{:});
			
			loc_ind_voxel = {i:i+1, j:j+1, k:k+1};
			
% 			% Take care of the objects touching boundaries.
% 			for v = Axis.elems
% 				if gt_array(v) == GT.dual && loc_ind_voxel{v}(1) == 1
% 					loc_ind_voxel{v} = loc_ind_voxel{v}(end);
% 				elseif gt_array(v) == GT.prim && loc_ind_voxel{v}(end) == grid3d.N(v)+1
% 					loc_ind_voxel{v} = loc_ind_voxel{v}(1);
% 				end
% 			end
			
			mat_ind_voxel = mat_ind_array(loc_ind_voxel{:});
			shape_ind_voxel = shape_ind_array(loc_ind_voxel{:});
			
			unique_mat_ind = unique(mat_ind_voxel(:));  % unique and sorted
			n_mat_ind = length(unique_mat_ind);
			assert(n_mat_ind >= 2);  % voxel is heterogeneous
			
			voxel = [l_voxel{Axis.x}(i:i+1); l_voxel{Axis.y}(j:j+1); l_voxel{Axis.z}(k:k+1)];  % voxel
			voxel_center = [l_voxel_center{Axis.x}(i), l_voxel_center{Axis.y}(j), l_voxel_center{Axis.z}(k)];  % voxel

			mat_bg_ind = unique_mat_ind(1);  % background material
			mat_bg = ind2mat_array(:, :, mat_bg_ind);
			for mat_ind_k = 2:n_mat_ind
				mat_fg_ind = unique_mat_ind(mat_ind_k);  % foreground material (which is added later)
				mat_fg = ind2mat_array(:, :, mat_fg_ind);

				where_mat_fg = (mat_ind_voxel == mat_fg_ind);
				shape_ind_fg = shape_ind_voxel(where_mat_fg);

				unique_shape_inds = unique(shape_ind_fg(:));  % column
				unique_shape_inds = unique_shape_inds.';  % row
				if length(unique_shape_inds) == 1
					unique_shape = ind2shape_array(unique_shape_inds);
				else
					assert(length(unique_shape_inds) >= 2);
					key = mat2str(unique_shape_inds);
					if union_shape_map.isKey(key)  % union shape already created
						unique_shape = union_shape_map(key);
					else % union shape not created yet
						unique_shape = UnionShape(ind2shape_array(unique_shape_inds));
						union_shape_map(key) = unique_shape;
					end
				end

% 				[rvol, ndir] = unique_shape.smoothing_params(voxel, voxel_center);
% 				[rvol, ndir] = unique_shape.smoothing_params(voxel);
				rvol = unique_shape.fill_factor(voxel);
				ndir = unique_shape.outward_normal(voxel);
				xyz0 = 'xyz0'; fprintf('%s%s at %s: rvol = %f, ndir = %s\n', char(ft), xyz0(w), mat2str(voxel_center), rvol, mat2str(ndir));
				
				mat_bg = smooth_tensor(mat_fg, mat_bg, rvol, ndir);
				if any(isnan(mat_bg(:)))  % possible if norm(ndir) = 0
					exception = MException('Maxwell:subpixel', 'Material parameter tensor has NaN.');
					throwAsCaller(exception);
				end
			end

			mat_smoothed = mat_bg;
			
			if w == 1+Axis.count
				ind_offdiag = (eye(Axis.count) == 0);
				if ft == FT.e
					eps_array(i, j, k, ind_offdiag) = mat_smoothed(ind_offdiag);
				else
					mu_array(i, j, k, ind_offdiag) = mat_smoothed(ind_offdiag);
				end
			else  % w == x, y, z
				if ft == FT.e
					eps_array(i, j, k, w, w) = mat_smoothed(w, w);
				else
					mu_array(i, j, k, w, w) = mat_smoothed(w, w);
				end
			end
		end
	end
end

eps_array = ipermute(eps_array, [3 4 5 1 2]);
mu_array = ipermute(mu_array, [3 4 5 1 2]);


function logical_array = find_hetero_voxels(grid3d, mat_ind_array)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" is instance of Grid3d.');
chkarg(istypesizeof(mat_ind_array, 'int', grid3d.N + 1), ...
	'"mat_ind_array" should be %d-by-%d-by-%d array with integer elements.', ...
	grid3d.N(Axis.x)+1, grid3d.N(Axis.y)+1, grid3d.N(Axis.z)+1);

logical_array = false(grid3d.N);
ind_shift = [0 1];

for dk = ind_shift
	for dj = ind_shift
		for di = ind_shift
			if isequal(di, dj, dk, 1)
				break
			end
			logical_array = logical_array ...
				| (mat_ind_array(2-di:end-di, 2-dj:end-dj, 2-dk:end-dk) ...
					~= mat_ind_array(1:end-1, 1:end-1, 1:end-1));
		end
	end
end

function mat_ind_array = attach_ghost_mat_ind(mat_ind_array, gt_array, bc_array)
% Augment the Fw array with ghost boundary elements in the v-axis.
% Surprisingly, this function combines attach_extra_E() and attach_extra_H().

chkarg(istypesizeof(mat_ind_array, 'int', [0 0 0]), '"gt_array" should be %dD array with integer elements.', Axis.count);
chkarg(istypesizeof(gt_array, 'GT', [1, Axis.count]), '"gt_array" should be length-%d row vector with GT as elements.', Axis.count);
chkarg(istypesizeof(bc_array, 'BC', [1, Axis.count]), '"bc_array" should be length-%d row vector with BC as elements.', Axis.count);

for w = Axis.elems
	ind = {':',':',':'};
	gt = gt_array(w);
	bc_w = bc_array(w);
	Nw = size(mat_ind_array, int(w));
	if gt == GT.prim
		if bc_w == BC.p
			ind{w} = 1;
		else  % bc_w == BC.e or BC.m
			ind{w} = Nw;
		end
		mat_ind_array = cat(int(w), mat_ind_array, mat_ind_array(ind{:}));
	else  % gt == GT.dual
		if bc_w == BC.p
			ind{w} = Nw;
		else  % bc_w == BC.e or BC.m
			ind{w} = 1;
		end
		mat_ind_array = cat(int(w), mat_ind_array(ind{:}), mat_ind_array);
	end
end

function [shape_ind_array, ind2shape_array] = attach_ghost_shape_ind(domain, shape_ind_array, ind2shape_array, gt_array, bc_array)
% Augment the Fw array with ghost boundary elements in the v-axis.
% Surprisingly, this function combines attach_extra_E() and attach_extra_H().

chkarg(istypesizeof(domain, 'Domain'), '"domain" shuold be instance of Domain.');
chkarg(istypesizeof(shape_ind_array, 'int', [0 0 0]), '"shape_array" should be %dD array with integer elements.', Axis.count);
chkarg(istypesizeof(ind2shape_array, 'Shape', [1 0]), '"ind2shape_array" should be row vector with Shape as elements.');
chkarg(istypesizeof(gt_array, 'GT', [1, Axis.count]), '"gt_array" should be length-%d row vector with GT as elements.', Axis.count);
chkarg(istypesizeof(bc_array, 'BC', [1, Axis.count]), '"bc_array" should be length-%d row vector with BC as elements.', Axis.count);

for w = Axis.elems
	ind = {':',':',':'};
	gt = gt_array(w);
	bc_w = bc_array(w);
	Nw = size(shape_ind_array, int(w));
		
	if gt == GT.prim
		sign = Sign.p;		
		if bc_w == BC.p
			ind{w} = 1;
		else  % bc_w == BC.e or BC.m
			ind{w} = Nw;
		end
	else  % gt == GT.dual
		sign = Sign.n;
		if bc_w == BC.p
			ind{w} = Nw;
		else  % bc_w == BC.e or BC.m
			ind{w} = 1;
		end
	end
	
	shape_ind_g_array = shape_ind_array(ind{:});  % ghost shapes
	unique_shape_inds = unique(shape_ind_g_array(:));  % column
	unique_shape_inds = unique_shape_inds.';  % row
	for shape_ind = unique_shape_inds
		shape = ind2shape_array(shape_ind);
		if isempty(shape.domain)
			shape.domain = domain;  % this adds domain to shape before ghost point because Shape is handle class
		end

		% Create a ghost shape.
		if bc_w == BC.p
			shape_g = shape.shift_domain(w, sign);
		else  % bc_w == BC.e or BC.m
			shape_g = shape.flip_domain(w, sign);
		end

		ind2shape_array = [ind2shape_array(1:end), shape_g];
		shape_ind_g = length(ind2shape_array);
		shape_ind_g_array(shape_ind_g_array == shape_ind) = shape_ind_g;
	end

	if gt == GT.prim
		shape_ind_array = cat(int(w), shape_ind_array, shape_ind_g_array);
	else  % gt == GT.dual
		shape_ind_array = cat(int(w), shape_ind_g_array, shape_ind_array);
	end
end
