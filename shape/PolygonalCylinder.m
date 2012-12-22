%% PolygonalCylinder
% Concrete subclass of <GenericCylinder.html |GenericCylinder|> representing a
% cylinder with a polygonal cross section.

%%% Description
% |PolygonalCylinder| is a Shape for a cylinder whose cross section is a
% polygon.  It is assumed that the polygon is convex and the mean of its
% vertices is inside the polygon.  The axis of the cylinder should be aligned
% with one of the axes of the Cartesian coordinate system.

%%% Construction
%  shape = PolygonalCylinder(normal_axis, height, normal_center, point_array)
%  shape = PolygonalCylinder(normal_axis, height, normal_center, point_array, dl_max)
% 
% *Input Arguments*
%
% * |normal_axis|: axis of the cylinder.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |height|: size of the cylinder along its axis.
% * |normal_center|: center of the axis of the cylinder.  It is a real number.
% * |radius|: radius of the circle
% * |point_array|: |n|-by-|2| array of the vertices of the polygon with |n|
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
% <CircularCylinder.html |CircularCylinder|>, <EllipticCylinder.html
% |EllipticCylinder|>, <SectoralCylinder.html |SectoralCylinder|>, <Shape.html
% |Shape|>, <maxwell_run.html |maxwell_run|>

classdef PolygonalCylinder < GenericCylinder

	methods
        function this = PolygonalCylinder(normal_axis, height, normal_center, point_array, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');
			chkarg(istypesizeof(normal_center, 'real'), '"normal_center" should be real.');
			chkarg(istypesizeof(point_array, 'real', [0 Dir.count]), ...
				'"point_array" should be %d-column matrix with real elements', Dir.count);

			Np = size(point_array, 1);
			p_array = point_array.';
			p_array = [p_array, p_array(:,1)];  % # of columns: Np+1
			c = mean(p_array, 2);
			pc_array = [p_array(1,:) - c(1); p_array(2,:) - c(2)];
			z_array = pc_array(1,:) + 1i * pc_array(2,:);
			theta_array = angle(z_array(2:end) ./ z_array(1:end-1));  % length: Np
			chkarg(all(theta_array > 0), ...
				'"point_array" should represent points arranged in counter-clockwise direction.');
			
			d_array = diff(p_array, [], 2);
			
			% Set up n_array, the array of outward unit vectors.
			n_array = [d_array(2,:); -d_array(1,:)];
			nn_array = sqrt(sum(n_array.^2));
			n_array = [n_array(1,:)./nn_array; n_array(2,:)./nn_array];  % # of columns: Np
			
			% Set up s_array, the array of distances from the center.
			s_array = sum(pc_array(:,1:end-1) .* n_array);  % length: Np
			chkarg(all(s_array > 0), 'mean point of vertices should be in polygon.');
			
			% lsf2d() can andle rho = [p q] with column vectors p and q.
			function level = lsf2d(rho)
				chkarg(istypesizeof(rho, 'real', [0, Dir.count]), ...
					'"rho" should be matrix with %d columns with real elements.', Dir.count);
				N = size(rho, 1);
				c_vec = repmat(c.', [N 1]);
				rc = rho - c_vec;
				rc = repmat(rc, [1 1 Np]);
				n_array_vec = reshape(n_array, [1 Dir.count Np]);
				n_array_vec = repmat(n_array_vec, [N 1]);
				s_array_vec = reshape(s_array, [1 1 Np]);
				s_array_vec = repmat(s_array_vec, [N 1]);
				level = 1 - max(sum(rc .* n_array_vec, 2) ./s_array_vec, [], 3);
			end

			[h, v, n] = cycle(normal_axis);
			lprim = cell(1, Axis.count);
			lprim{h} = p_array(1, 1:end-1);
			lprim{v} = p_array(2, 1:end-1);
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

