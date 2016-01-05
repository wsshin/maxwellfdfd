function [eps_array, mu_array] = mean_subpixel(grid3d, ge, domain, eps_imat_cell, mu_imat_cell, eps_ishape_cell, mu_ishape_cell, ind2eps_array, ind2mu_array, ind2shape_array, eps_array, mu_array)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');
chkarg(istypesizeof(domain, 'Domain'), '"domain" shuold be instance of Domain.');
chkarg(istypesizeof(eps_imat_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"eps_imat_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(mu_imat_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"mu_imat_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(eps_ishape_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"eps_ishape_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});
chkarg(istypesizeof(mu_ishape_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"mu_ishape_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer as elements.', ...
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

% Find heterogeneous voxels.
[I, J, K] = ndgrid(1:grid3d.N(Axis.x), 1:grid3d.N(Axis.y), 1:grid3d.N(Axis.z));
ijk_hetero_cell = cell(FT.count, Axis.count+1);
voxel_hetero_cell = cell(FT.count, Axis.count+1);
imat_hetero_cell = cell(FT.count, Axis.count+1);
ishape_hetero_cell = cell(FT.count, Axis.count+1);
for ft = FT.elems  % eps or mu
	for w = [1+Axis.count int(Axis.elems)]
		if ft == FT.e  % eps
			gt = alter(ge);  % for smoothing eps, voxels to examine have H-fields at vertices
			imat_array = eps_imat_cell{w};
			ishape_array = eps_ishape_cell{w};
		else  % mu
			gt = ge;  % for smoothing mu, voxels to examine have E-fields at vertices
			imat_array = mu_imat_cell{w};
			ishape_array = mu_ishape_cell{w};
		end

		gt_array = gt(ones(1, Axis.count));
		if w ~= 1+Axis.count  % w == x, y, z
			gt_array(w) = alter(gt);
		end

		l_voxel = grid3d.lg(Axis.elems + Axis.count*subsindex(gt_array));  % locations of vertices of voxels
		
		imat_array = attach_ghost_mat_ind(imat_array, gt_array, grid3d.bc);  % (Nx+1) x (Ny+1) x (Nz+1)
		[ishape_array, ind2shape_array] = attach_ghost_shape_ind(domain, ishape_array, ind2shape_array, gt_array, grid3d.bc);  % (Nx+1) x (Ny+1) x (Nz+1)

		is_hetero_array = find_hetero_voxels(grid3d, imat_array);  % imagesc(is_hetero_mat_array) shows locations to smooth

		ijk_hetero = [I(is_hetero_array), J(is_hetero_array), K(is_hetero_array)];  % n_hetero x 3
		
		n_hetero = size(ijk_hetero, 1);  % number of hetero voxels
		voxel_hetero_array = NaN(Axis.count, Sign.count, n_hetero);	
		for s = Sign.elems
			i_shape = ijk_hetero(:, Axis.x).' + subsindex(s);  % row
			js = ijk_hetero(:, Axis.y).' + subsindex(s);  % row
			ks = ijk_hetero(:, Axis.z).' + subsindex(s);  % row
					
			voxel_hetero_array(:, s, :) = [l_voxel{Axis.x}(i_shape); l_voxel{Axis.y}(js); l_voxel{Axis.z}(ks)];  % Axis.count x n_hetero
		end
		
		imat_hetero_array = NaN(Sign.count, Sign.count, Sign.count, n_hetero);
		ishape_hetero_array = NaN(Sign.count, Sign.count, Sign.count, n_hetero);	
		for sx = Sign.elems
			for sy = Sign.elems
				for sz = Sign.elems
					i_shape = ijk_hetero(:, Axis.x).' + subsindex(sx);  % row
					js = ijk_hetero(:, Axis.y).' + subsindex(sy);  % row
					ks = ijk_hetero(:, Axis.z).' + subsindex(sz);  % row
					
					imat_hetero_array(sx, sy, sz, :) = imat_array(sub2ind(grid3d.Ng, i_shape, js, ks));  % n_hetero x Axis.count
					ishape_hetero_array(sx, sy, sz, :) = ishape_array(sub2ind(grid3d.Ng, i_shape, js, ks));  % n_hetero x Axis.count
				end
			end
		end
				
		ijk_hetero_cell{ft, w} = ijk_hetero;
		voxel_hetero_cell{ft, w} = voxel_hetero_array;
		imat_hetero_cell{ft, w} = imat_hetero_array;
		ishape_hetero_cell{ft, w} = ishape_hetero_array;
	end
end

% Replace multiple shapes corresponding to the same material in a voxel with a
% union shape.
ishape_union_map = containers.Map();
for ft = FT.elems  % eps or mu
	for w = [1+Axis.count int(Axis.elems)]
		ijk_hetero = ijk_hetero_cell{ft, w};
		imat_hetero_array = imat_hetero_cell{ft, w};
		ishape_hetero_array = ishape_hetero_cell{ft, w};
		
		n_hetero = size(ijk_hetero, 1);  % number of hetero voxels
		for i_hetero = 1:n_hetero
			imat_voxel = imat_hetero_array(:,:,:,i_hetero);  % materials inside voxel
			ishape_voxel = ishape_hetero_array(:,:,:,i_hetero);  % shapes inside voxel
			
			unique_imat = unique(imat_voxel(:));  % unique and sorted
			n_imat = length(unique_imat);
			assert(n_imat >= 2);  % voxel is heterogeneous

			for k_imat = 2:n_imat  % 1st material is background material, so ignored
				imat_fg = unique_imat(k_imat);  % foreground material (which is added later)

				where_imat_fg = (imat_voxel == imat_fg);
				ishape_fg = ishape_voxel(where_imat_fg);

				unique_ishapes = unique(ishape_fg(:));  % column and sorted
				unique_ishapes = unique_ishapes.';  % row and sorted
				if length(unique_ishapes) >= 2
					key = mat2str(unique_ishapes);  % generate string as key
					if ishape_union_map.isKey(key)  % union shape already created
						ishape_union = ishape_union_map(key);
					else % union shape yet to be created
						shape_union = UnionShape(ind2shape_array(unique_ishapes));
						ind2shape_array = [ind2shape_array(1:end), shape_union];
						ishape_union = length(ind2shape_array);
						ishape_union_map(key) = ishape_union;
					end
					
					ishape_voxel(where_imat_fg) = ishape_union;
					ishape_hetero_array(:,:,:,i_hetero) = ishape_voxel;
				end
			end

		end
		
		ishape_hetero_cell{ft, w} = ishape_hetero_array;
	end
end

% Construct a cell array mapping shape indices to the voxels intersecting the shapes.
n_shape = length(ind2shape_array);
ishape2voxel_cell = cell(1, n_shape);
for ft = FT.elems  % eps or mu
	for w = [1+Axis.count int(Axis.elems)]
		ijk_hetero = ijk_hetero_cell{ft, w};
		voxel_hetero_array = voxel_hetero_cell{ft, w};
		ishape_hetero_array = ishape_hetero_cell{ft, w};

		n_hetero = size(ijk_hetero, 1);  % number of hetero voxels
		ishape_hetero_array = reshape(ishape_hetero_array, Sign.count^Axis.count, n_hetero);
		
		for i_shape = 1:n_shape
			ishape2voxel_array = ishape2voxel_cell{i_shape};
			
			has_voxel_shape = any(ishape_hetero_array == i_shape);  % row
			voxels_with_shape = voxel_hetero_array(:, :, has_voxel_shape);
			ishape2voxel_array = cat(3, ishape2voxel_array, voxels_with_shape);
			
			ishape2voxel_cell{i_shape} = ishape2voxel_array;
		end
	end
end

% Calculate the fill factors and outward normals of the shapes in hetero voxels.
ishape2rvol_cell = cell(1, n_shape);
ishape2ndir_cell = cell(1, n_shape);
for i_shape = 1:n_shape
	shape = ind2shape_array(i_shape);
	voxel_array = ishape2voxel_cell{i_shape};  % scatter(voxel_array(Axis.x,:), voxel_array(Axis.y,:)) shows voxel distribution
	if ~isempty(voxel_array)
		tic; ishape2rvol_cell{i_shape} = shape.fill_factor(voxel_array); toc;
		tic; ishape2ndir_cell{i_shape} = shape.outward_normal(voxel_array); toc;
	end
end

% Perform subpixel smoothing using the precalculated fill factors and outward
% normals.
shape_voxel_counter = zeros(1, n_shape);
for ft = FT.elems  % eps or mu
	if ft == FT.e  % eps
		ind2mat_array = ind2eps_array;
	else  % mu
		ind2mat_array = ind2mu_array;
	end

	for w = [1+Axis.count int(Axis.elems)]
		ijk_hetero = ijk_hetero_cell{ft, w};
		imat_hetero_array = imat_hetero_cell{ft, w};
		ishape_hetero_array = ishape_hetero_cell{ft, w};
		
		n_hetero = size(ijk_hetero, 1);  % number of hetero voxels
		for i_hetero = 1:n_hetero
			imat_voxel = imat_hetero_array(:,:,:,i_hetero);  % materials inside voxel
			ishape_voxel = ishape_hetero_array(:,:,:,i_hetero);  % shapes inside voxel
			
			unique_imat = unique(imat_voxel(:));  % unique and sorted
			n_imat = length(unique_imat);
			assert(n_imat >= 2);  % voxel is heterogeneous

			imat_bg = unique_imat(1);  % background material
			mat_bg = ind2mat_array(:, :, imat_bg);
			for k_imat = 2:n_imat  % 1st material is background material, so ignored
				imat_fg = unique_imat(k_imat);  % foreground material (which is added later)
				mat_fg = ind2mat_array(:, :, imat_fg);

				where_imat_fg = (imat_voxel == imat_fg);
				ishape_fg = ishape_voxel(where_imat_fg);

				i_shape = unique(ishape_fg(:));  % column and sorted
				assert(length(i_shape) == 1);  % union shapes already generated
				
				c_shape = shape_voxel_counter(i_shape);  % counter
				c_shape = c_shape + 1;
				shape_voxel_counter(i_shape) = c_shape;
				
				rvol = ishape2rvol_cell{i_shape}(c_shape);
				ndir = ishape2ndir_cell{i_shape}(:, c_shape);
				
% 				voxel = ishape2voxel_cell{i_shape}(:,:,c_shape);
% 				voxel_center = mean(voxel.');
% 				xyz0 = 'xyz0'; fprintf('%s%s at %s: rvol = %f, ndir = %s\n', char(ft), xyz0(w), mat2str(voxel_center), rvol, mat2str(ndir));

				mat_bg = smooth_tensor(mat_fg, mat_bg, rvol, ndir);
				if any(isnan(mat_bg(:)))  % possible if norm(ndir) = 0
					exception = MException('Maxwell:subpixel', 'Material parameter tensor has NaN.');
					throwAsCaller(exception);
				end
			end
			mat_smoothed = mat_bg;

			ijk = num2cell(ijk_hetero(i_hetero, :));
			if w == 1+Axis.count
				ind_offdiag = (eye(Axis.count) == 0);
				if ft == FT.e
					eps_array(ijk{:}, ind_offdiag) = mat_smoothed(ind_offdiag);
				else
					mu_array(ijk{:}, ind_offdiag) = mat_smoothed(ind_offdiag);
				end
			else  % w == x, y, z
				if ft == FT.e
					eps_array(ijk{:}, w, w) = mat_smoothed(w, w);
				else
					mu_array(ijk{:}, w, w) = mat_smoothed(w, w);
				end
			end
		end
	end
end

% for i_shape = 1:n_shape
% 	assert(shape_voxel_counter{i_shape} == size(ishape2voxel_cell{i_shape}, 3));
% end

eps_array = ipermute(eps_array, [3 4 5 1 2]);
mu_array = ipermute(mu_array, [3 4 5 1 2]);

function logical_array = find_hetero_voxels(grid3d, imat_array)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" is instance of Grid3d.');
chkarg(istypesizeof(imat_array, 'int', grid3d.Ng), ...
	'"imat_array" should be %d-by-%d-by-%d array with integer elements.', grid3d.Ngcell{:});

logical_array = false(grid3d.N);
ind_shift = [0 1];

for dk = ind_shift
	for dj = ind_shift
		for di = ind_shift
			if isequal(di, dj, dk, 1)
				break
			end
			logical_array = logical_array ...
				| (imat_array(2-di:end-di, 2-dj:end-dj, 2-dk:end-dk) ...
					~= imat_array(1:end-1, 1:end-1, 1:end-1));
		end
	end
end

function imat_array = attach_ghost_mat_ind(imat_array, gt_array, bc_array)
% Augment the Fw array with ghost boundary elements in the v-axis.
% Surprisingly, this function combines attach_extra_E() and attach_extra_H().

chkarg(istypesizeof(imat_array, 'int', [0 0 0]), '"imat_array" should be %dD array with integer elements.', Axis.count);
chkarg(istypesizeof(gt_array, 'GT', [1, Axis.count]), '"gt_array" should be length-%d row vector with GT as elements.', Axis.count);
chkarg(istypesizeof(bc_array, 'BC', [1, Axis.count]), '"bc_array" should be length-%d row vector with BC as elements.', Axis.count);

for w = Axis.elems
	ind = {':',':',':'};
	gt = gt_array(w);
	bc_w = bc_array(w);
	Nw = size(imat_array, int(w));
	if gt == GT.prim
		if bc_w == BC.p
			ind{w} = 1;
		else  % bc_w == BC.e or BC.m
			ind{w} = Nw;
		end
		imat_array = cat(int(w), imat_array, imat_array(ind{:}));
	else  % gt == GT.dual
		if bc_w == BC.p
			ind{w} = Nw;
		else  % bc_w == BC.e or BC.m
			ind{w} = 1;
		end
		imat_array = cat(int(w), imat_array(ind{:}), imat_array);
	end
end

function [ishape_array, ind2shape_array] = attach_ghost_shape_ind(domain, ishape_array, ind2shape_array, gt_array, bc_array)
% Augment the Fw array with ghost boundary elements in the v-axis.
% Surprisingly, this function combines attach_extra_E() and attach_extra_H().

chkarg(istypesizeof(domain, 'Domain'), '"domain" shuold be instance of Domain.');
chkarg(istypesizeof(ishape_array, 'int', [0 0 0]), '"ishape_array" should be %dD array with integer elements.', Axis.count);
chkarg(istypesizeof(ind2shape_array, 'Shape', [1 0]), '"ind2shape_array" should be row vector with Shape as elements.');
chkarg(istypesizeof(gt_array, 'GT', [1, Axis.count]), '"gt_array" should be length-%d row vector with GT as elements.', Axis.count);
chkarg(istypesizeof(bc_array, 'BC', [1, Axis.count]), '"bc_array" should be length-%d row vector with BC as elements.', Axis.count);

for w = Axis.elems
	ind = {':',':',':'};
	gt = gt_array(w);
	bc_w = bc_array(w);
	Nw = size(ishape_array, int(w));

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

	ishape_g_array = ishape_array(ind{:});  % ghost shapes
	unique_ishapes = unique(ishape_g_array(:));  % column
	unique_ishapes = unique_ishapes.';  % row
	for i_shape = unique_ishapes
		shape = ind2shape_array(i_shape);
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
		ishape_g = length(ind2shape_array);
		ishape_g_array(ishape_g_array == i_shape) = ishape_g;
	end

	if gt == GT.prim
		ishape_array = cat(int(w), ishape_array, ishape_g_array);
	else  % gt == GT.dual
		ishape_array = cat(int(w), ishape_g_array, ishape_array);
	end
end
