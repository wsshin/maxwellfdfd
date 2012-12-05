%% periodize_shape
% Constructs a periodic array of a given shape.

%%% Syntax
%  shape_array = periodize_shape(shape, primitive_vec_cell, boundshape)

%%% Parameters
% *Input*
%
% * |shape|: instance of <Shape.html |Shape|> to be repeated.
% * |primitive_vec_cell|: three primitive lattice vectors in the format of |{[x1
% y1 z1], [x2 y2 z2], [x3 y3 z3]}|.  
% * |boundshape|: shape that contains the centers of the translated shapes.  It
% is used to control the overall extent of the periodic array.
%
% *Output*
%
% * |shape_array|: array of instances of <Shape.html |Shape|>.

%%% Example
%   % Create a shape for a slab.
%   t = 500;  % thickness
%   slab = Box([-5000 5000; -3000 3000; 0 t]);
%
%   % Create an array of circular cylinders.
%   r = 100;  % radius
%   rod = CircularCylinder(Axis.z, t, [0 0 t/2], r);
%   rod_array = periodize_shape(rod, {[400 0 0], [0 400 0], [0 0 t]}, slab)
%   
%   % Use the constructed shape in maxwell_run().
%   gray = [0.5 0.5 0.5];
%   [E, H] = maxwell_run(...
%               {INITIAL ARGUMENTS}, ...
%               'OBJ', ...
%                   {'silicon', gray, 11.56}, slab, ...
%                   {'vacuum', 'none', 1.0}, rod_array, ...
%               {REMAINING ARGUMENTS});

%%% See Also
% <Shape.html |Shape|>, <maxwell_run.html |maxwell_run|>

function shape_array = periodize_shape(shape, primitive_vec_cell, boundshape)

chkarg(istypesizeof(shape, 'Shape'), '"shape" should be instance of Shape.');
chkarg(istypesizeof(primitive_vec_cell, 'realcell', [1 Axis.count], [1 Axis.count]), ...
	'"primitive_vec_cell" should be length-%d row cell array whose each element is length-%d row vector with real elements.', Axis.count, Axis.count);
chkarg(istypesizeof(boundshape, 'Shape'), '"boundshape" should be instance of Shape.')

a = primitive_vec_cell;  % "a" is conventional symbol for primitive vectors
shape0 = shape;

lsf0 = shape0.lsf;  % level set function
lprim0 = shape0.lprim;  % primary grid planes
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

nshape = prod((nmax-nmin)+1);  % tentative number of shapes
% for nz = nmin(Axis.z):nmax(Axis.z)
% 	for ny = nmin(Axis.y):nmax(Axis.y)
% 		for nx = nmin(Axis.x):nmax(Axis.x)
% 			n = [nx ny nz];
% 			R = n * A;  % lattice vector
% 			if boundshape.contains(c0 + R)
% 				nshape = nshape + 1;
% 			end
% 		end
% 	end
% end
shape_array = Shape.empty(0, nshape);

ishape = 0;
for nz = nmin(Axis.z):nmax(Axis.z)
	for ny = nmin(Axis.y):nmax(Axis.y)
		for nx = nmin(Axis.x):nmax(Axis.x)
			n = [nx ny nz];
			R = n * A;  % lattice vector
			if boundshape.contains(c0 + R)
				lprim = cell(1, Axis.count);
				for w = Axis.elems
					lprim{w} = lprim0{w} + R(w);
				end
				lsf = @(r) lsf0(r - repmat(R, [size(r,1), 1]));
				shape = Shape(lprim, lsf, dl_max0);
				ishape = ishape + 1;
				shape_array(ishape) = shape;
			end
		end
	end
end
