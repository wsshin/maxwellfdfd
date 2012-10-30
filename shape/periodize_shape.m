function shape_array = periodize_shape(shape, primitive_vec_cell, boundshape)

chkarg(istypesizeof(shape, 'Shape'), '"shape" should be instance of Shape.');
chkarg(istypesizeof(primitive_vec_cell, 'realcell', [1 Axis.count], [1 Axis.count]), ...
	'"primitive_vec_cell" should be length-%d row cell array whose each element is length-%d row vector with real elements.', Axis.count, Axis.count);
chkarg(istypesizeof(boundshape, 'Shape'), '"boundshape" should be instance of Shape.')

a = primitive_vec_cell;  % "a" is conventional symbol for primitive vectors
shape0 = shape;

lsf0 = shape0.lsf;  % level set function
box0 = shape0.bound;  % circumbox
c0 = shape0.cb_center;  % center of circumbox
dl_max0 = shape0.dl_max;  % dl_max

box_bound = boundshape.bound;  % circumbox of boundshape
n_corners = Sign.count^Axis.count;  % 8
corners = NaN(n_corners, Axis.count);  % 8-by-3
ic = 1;
for sz = Sign.elems
	for sy = Sign.elems
		for sx = Sign.elems
			corners(ic, Axis.x) = box_bound(Axis.x, sx);
			corners(ic, Axis.y) = box_bound(Axis.y, sy);
			corners(ic, Axis.z) = box_bound(Axis.z, sz);
			ic = ic + 1;
		end
	end
end


A = [a{1}; a{2}; a{3}];

% Calculate all n = [nx ny nz] that make n*A + c0 corners.  These are extreme n's.
n_ext = (corners - repmat(c0, [n_corners, 1])) / A;  % n * A + c0 = corner

% nmax = round(max(n_ext-eps, [], 1));  % round() instead of floor() to account for round-off error; -eps is to round down +X.5
% nmin = round(min(n_ext+eps, [], 1));  % round() instead of ceil() to account for round-off error; +eps is to round up -X.5
nmax = floor(max(n_ext-eps, [], 1));
nmin = ceil(min(n_ext+eps, [], 1));

shape_array = [];
for nz = nmin(Axis.z):nmax(Axis.z)
	for ny = nmin(Axis.y):nmax(Axis.y)
		for nx = nmin(Axis.x):nmax(Axis.x)
			n = [nx ny nz];
			R = n * A;  % lattice vector
			if boundshape.contains(c0 + R)
				circumbox = box0 + [R; R].';
				lsf = @(r) lsf0(r - repmat(R, [size(r,1), 1]));
				if isnan(dl_max0)
					shape = Shape(circumbox, lsf);
				else
					shape = Shape(circumbox, lsf, dl_max0);
				end
				shape_array = [shape_array(1:end), shape];
			end
		end
	end
end
