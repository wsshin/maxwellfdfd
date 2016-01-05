function [eps_imat_cell, mu_imat_cell, eps_ishape_cell, mu_ishape_cell, ind2eps_array, ind2mu_array, ind2shape_array, eps_array, mu_array] = ...
	assign_mat_ind(grid3d, ge, obj_array, eps_imat_cell, mu_imat_cell, eps_ishape_cell, mu_ishape_cell, ind2eps_array, ind2mu_array, ind2shape_array, eps_array, mu_array)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');
chkarg(istypesizeof(obj_array, 'Object', [1 0]), ...
	'"obj_array" should be row vector of instances of Object.');

if nargin < 4  % no eps_ind_cell
	eps_imat_cell = {NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N)};
end
chkarg(istypesizeof(eps_imat_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"eps_imat_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});

if nargin < 5  % no mu_ind_cell
	mu_imat_cell = {NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N)};
end
chkarg(istypesizeof(mu_imat_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"mu_imat_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});

if nargin < 6  % no eps_shape_ind_cell
	eps_ishape_cell = {NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N)};
end
chkarg(istypesizeof(eps_ishape_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"eps_ishape_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});

if nargin < 7  % no mu_shape_ind_cell
	mu_ishape_cell = {NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N), NaN(grid3d.N)};
end
chkarg(istypesizeof(mu_ishape_cell, 'intcell', [1, 1+Axis.count], grid3d.N), ...
	'"mu_ishape_cell" should be length-%d cell array whose each element is %d-by-%d-by-%d array with integer elements.', ...
	Axis.count, grid3d.Ncell{:});

if nargin < 8  % no ind2eps_array
	ind2eps_array = obj_array(1).material.eps;  % obj_array always has element
end
chkarg(istypesizeof(ind2eps_array, 'complex', [Axis.count Axis.count 0]), ...
	'"ind2eps_array" should be %d-by-%d-by-n array with complex elements.', Axis.count, Axis.count);

if nargin < 9  % no ind2mu_array
	ind2mu_array = obj_array(1).material.mu;  % obj_array always has element
end
chkarg(istypesizeof(ind2mu_array, 'complex', [Axis.count Axis.count 0]), ...
	'"ind2mu_array" should be %d-by-%d-by-n array with complex elements.', Axis.count, Axis.count);

if nargin < 10  % no ind2shape_array
	ind2shape_array = Shape.empty();
end
chkarg(isempty(ind2shape_array) || istypesizeof(ind2shape_array, 'Shape', [1 0]), ...
	'"ind2shape_array" should be empty, or row vector with Shape as elements.');

if nargin < 11  % no eps_array
	eps_array = NaN(Axis.count, Axis.count, grid3d.Ncell{:});
end
chkarg(istypesizeof(eps_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
	'"eps_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});

if nargin < 12  % no mu_array
	mu_array = NaN(Axis.count, Axis.count, grid3d.Ncell{:});
end
chkarg(istypesizeof(mu_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
	'"mu_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});

ind_bound = cell(1, Axis.count);  % indices
i_eps = size(ind2eps_array, 3);
i_mu = size(ind2mu_array, 3);
i_shape = numel(ind2shape_array);
for obj = obj_array
	shape = obj.shape;
	material = obj.material;

	% Below, even if the current material has been used before, treat it as a
	% different material and assign a different material index to it if other
	% material was used between these two same materials.  This mimicks the way
	% we use materials in 'OBJ' parameter group in maxwell_run().
	if ~isequal(ind2eps_array(:,:,end), material.eps)
		ind2eps_array = cat(3, ind2eps_array, material.eps);
		i_eps = i_eps + 1;  % i_eps == size(ind2eps_cell, 3)
	end
	
	if ~isequal(ind2mu_array(:,:,end), material.mu)
		ind2mu_array = cat(3, ind2mu_array, material.mu);
		i_mu = i_mu + 1;  % i_mu = size(ind2mu_cell, 3);
	end
	
	ind2shape_array = [ind2shape_array(1:end), shape];
	i_shape = i_shape + 1;
	
	for ft = FT.elems  % E or H
		if ft == FT.e
			gt = ge;

		else  % ft == FT.h
			gt = alter(ge);
		end
		
		for w = [1+Axis.count, int(Axis.elems)]  % assign material index to Ew or Hw locations
			gt_array = gt(ones(1, Axis.count));  % grid type for Fw locations
			if w ~= 1+Axis.count  % w == x, y, z
				gt_array(w) = alter(gt);
			end

			l_Fw = grid3d.l(Axis.elems + Axis.count*subsindex(gt_array));  % locations for Fw

			for v = Axis.elems
				bn = shape.bound(v, Sign.n);
				bp = shape.bound(v, Sign.p);
				in = find(l_Fw{v} >= bn, 1, 'first');
				ip = find(l_Fw{v} <= bp, 1, 'last');
				ind_bound{v} = in:ip;
			end
			
			if istypesizeof(shape, 'Box')
				if ft == FT.e
					if w == 1+Axis.count  % set off-diagonal entries of eps
						for r = Axis.elems
							c = cycle(r);
							while c ~= r
								eps_array(r, c, ind_bound{:}) = material.eps(r, c);
								c = cycle(c);
							end
						end
					else  % w == x, y, z; set diagonal entries of eps
						eps_array(w, w, ind_bound{:}) = material.eps(w, w);
					end
					
					% Ew locations are examined for smoothing mu.
					mu_imat_cell{w}(ind_bound{:}) = i_mu;
					mu_ishape_cell{w}(ind_bound{:}) = i_shape;
				else  % ft == FT.h
					if w == 1+Axis.count  % set off-diagonal entries of mu
						for r = Axis.elems
							c = cycle(r);
							while c ~= r
								mu_array(r, c, ind_bound{:}) = material.mu(r, c);
								c = cycle(c);
							end
						end
					else  % w == x, y, z; set diagonal entries of mu
						mu_array(w, w, ind_bound{:}) = material.mu(w, w);
					end
					
					% Hw locations are examined for smoothing eps.
					eps_imat_cell{w}(ind_bound{:}) = i_eps;
					eps_ishape_cell{w}(ind_bound{:}) = i_shape;
				end
			else  % shape is not a Box
				[X, Y, Z] = ndgrid(l_Fw{Axis.x}(ind_bound{Axis.x}), l_Fw{Axis.y}(ind_bound{Axis.y}), l_Fw{Axis.z}(ind_bound{Axis.z}));
				is_in = shape.contains(X, Y, Z);
				
				ind_logical = false(grid3d.N);  % logical indices
				ind_logical(ind_bound{:}) = is_in;
				
				if ft == FT.e
					if w == 1+Axis.count  % set off-diagonal entries of eps
						for r = Axis.elems
							c = cycle(r);
							while c ~= r
								eps_array(r, c, ind_logical) = material.eps(r, c);
								c = cycle(c);
							end
						end
					else  % w == x, y, z; set diagonal entries of eps
						eps_array(w, w, ind_logical) = material.eps(w, w);
					end
					
					% Ew locations are examined for smoothing mu.
					mu_imat_cell{w}(ind_logical) = i_mu;
					mu_ishape_cell{w}(ind_logical) = i_shape;
				else  % ft == FT.h
					if w == 1+Axis.count  % set off-diagonal entries of mu
						for r = Axis.elems
							c = cycle(r);
							while c ~= r
								mu_array(r, c, ind_logical) = material.mu(r, c);
								c = cycle(c);
							end
						end
					else  % w == x, y, z; set diagonal entries of mu
						mu_array(w, w, ind_logical) = material.mu(w, w);
					end
					
					% Hw locations are examined for smoothing eps.
					eps_imat_cell{w}(ind_logical) = i_eps;
					eps_ishape_cell{w}(ind_logical) = i_shape;
				end				
			end
		end
	end
end

% Check if eps_shape_cell{w} and mu_shape_cell{w} have the correct size; they
% were initially empty arrays.
for w = Axis.elems
	assert((grid3d.N(Axis.z)==1 && isequal(size(eps_ishape_cell{w}), grid3d.N([Axis.x Axis.y]))) ...
		|| isequal(size(eps_ishape_cell{w}), grid3d.N));
	assert(isequal(size(eps_ishape_cell{w}), size(mu_ishape_cell{w})));
end
