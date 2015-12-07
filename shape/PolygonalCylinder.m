%% PolygonalCylinder
% Concrete subclass of <GenericCylinder.html |GenericCylinder|> representing a
% cylinder with a polygonal cross section.

%%% Description
% |PolygonalCylinder| is a Shape for a cylinder whose cross section is a
% polygon.  It is assumed that the polygon is convex and the mean of its
% vertices is inside the polygon.  The axis of the cylinder should be aligned
% with one of the axes of the Cartesian coordinate system.

%%% Construction
%  shape = PolygonalCylinder(normal_axis, height, normal_center, vertex_array)
%  shape = PolygonalCylinder(normal_axis, height, normal_center, vertex_array, dl_max)
% 
% *Input Arguments*
%
% * |normal_axis|: axis of the cylinder.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |height|: size of the cylinder along its axis.
% * |normal_center|: center of the axis of the cylinder.  It is a real number.
% * |radius|: radius of the circle
% * |vertex_array|: |n|-by-|2| array of the vertices of the polygon with |n|
% vertices.  For |normal_axis == Axis.y|, it is in the format of |[z1 x1; ...;
% zn xn]|.
% * |dl_max|: maximum grid size allowed in the cylinder.  It can be either |[dx
% dy dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned,
% |dl_max = Inf| is used.

%%% Example
%   % Create an instance of SectoralCylinder.
%   shape = PolygonalCylinder(Axis.z, 100, 50, [0 0; 50 0; 25 25*sqrt(3)]);  % cross section is regular triangle
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <CircularCylinder.html |CircularCylinder|>, <CircularShellCylinder.html
% |CircularShellCylinder|>, <EllipticCylinder.html
% |EllipticCylinder|>, <SectoralCylinder.html |SectoralCylinder|>, <Shape.html
% |Shape|>, <maxwell_run.html |maxwell_run|>

classdef PolygonalCylinder < GenericCylinder

	methods
        function this = PolygonalCylinder(normal_axis, height, normal_center, vertex_array, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');
			chkarg(istypesizeof(normal_center, 'real'), '"normal_center" should be real.');
			chkarg(istypesizeof(vertex_array, 'real', [0 Dir.count]), ...
				'"vertex_array" should be %d-column matrix with real elements', Dir.count);

			v_array = vertex_array;  % shorthand notation
			
			Nv = size(v_array, 1);  % number of vertices
			c = mean(v_array);  % center of mass
			
			vc_array = [v_array; v_array(1,:)];  % # of vertices: Nv+1
			vc_array = bsxfun(@minus, vc_array, c);  % vertices in center-of-mass coordinates
			z_array = vc_array(:,1) + 1i * vc_array(:,2);  % coordinates in complex plane
			theta_array = angle(z_array(2:end) ./ z_array(1:end-1));  % length: Nv
			chkarg(all(theta_array > 0), ...
				'"vertex_array" should represent points arranged in counter-clockwise direction.');
			
			s_array = diff(vc_array);  % array of side vectors
			
			% Set up n_array, the array of outward unit vectors.
			n_array = [s_array(:,2), -s_array(:,1)];  % unnormalized outward normal directions
			ln_array = sqrt(sum(n_array.^2, 2));  % lengths of normal vectors
			n_array = bsxfun(@rdivide, n_array, ln_array);  % # of normal directions: Nv
			
			% Set up d_array, the array of distances from the center.
			d_array = sum(vc_array(1:end-1,:) .* n_array, 2);  % projections to normal directions (# of distances: Nv)
			chkarg(all(d_array > 0), 'mean point of vertices should be in polygon.');
			
			% lsf2d() can andle rho = [p q] with column vectors p and q.
			function level = lsf2d(p, q)
				chkarg(istypeof(p, 'real'), '"p" should be array with real elements.');
				chkarg(istypeof(q, 'real'), '"q" should be array with real elements.');
				chkarg(isequal(size(p), size(q)), '"p" and "q" should have same size.');

				lc = {p - c(Dir.h), q - c(Dir.v)};  % locations in center-of-mass coordinates
				
				level = -Inf(size(p));
				for i = 1:Nv
					leveli = zeros(size(p));
					for d = Dir.elems
						leveli = leveli + lc{d} .* n_array(i,d);
					end
					leveli = leveli ./ d_array(i);
					level = max(level, leveli);
				end
				level = 1 - level;

% 				N = size(rho, 1);
% 				c_vec = repmat(c.', [N 1]);
% 				rc = rho - c_vec;
% 				rc = repmat(rc, [1 1 Nv]);
% 				n_array_vec = reshape(n_array, [1 Dir.count Nv]);
% 				n_array_vec = repmat(n_array_vec, [N 1]);
% 				d_array_vec = reshape(d_array, [1 1 Nv]);
% 				d_array_vec = repmat(d_array_vec, [N 1]);
% 				level = 1 - max(sum(rc .* n_array_vec, 2) ./d_array_vec, [], 3);
			end

			[h, v, n] = cycle(normal_axis);
			lprim = cell(1, Axis.count);
			lprim{h} = v_array(:, Dir.h).';
			lprim{v} = v_array(:, Dir.v).';
			lprim{n} = [-height height]/2 + normal_center;
			
			if nargin < 5  % no dl_max
				super_args = {normal_axis, @lsf2d, lprim};
			else
				super_args = {normal_axis, @lsf2d, lprim, dl_max};
			end
			
			this = this@GenericCylinder(super_args{:});
		end
	end
end

