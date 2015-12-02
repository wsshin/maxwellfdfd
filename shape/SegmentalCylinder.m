%% SegmentalCylinder
% Concrete subclass of <GenericCylinder.html |GenericCylinder|>
% representing a cylinder with a cross section of a circular segment.

%%% Description
% |SegmentalCylinder| represents the shape of a segmental cylinder.  Its
% cross section is a circular segment (region surrounded by an arc and the
% chord connecting its end points).  The axis of the cylinder should be
% aligned with one of the axes of the Cartesian coordinate system.

%%% Construction
%  shape = SegmentalCylinder(normal_axis, height, center, radius, theta, d_theta)
%  shape = SegmentalCylinder(normal_axis, height, center, radius, theta, d_theta, dl_max)
% 
% *Input Arguments*
%
% * |normal_axis|: axis of the cylinder.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |height|: size of the cylinder along its axis.
% * |center|: center of the cylinder in the format of |[x y z]|.  For
% |normal_axis = Axis.z|, |(x, y)| is the coordinate of the center of the
% circle.
% * |radius|: radius of the circle
% * |theta|: beginning angle of the segment in radian
% * |d_theta|: angular width of the segment in radian between -|pi| and |pi|
% (the angle of a segment should not be reflex)
% * |dl_max|: maximum grid size allowed in the cylinder.  It can be either |[dx
% dy dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned,
% |dl_max = Inf| is used.

%%% Example
%   % Create an instance of SegmentalCylinder.
%   shape = SegmentalCylinder(Axis.z, 100, [0 0 50], 50, pi/6, pi/3);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <SectoralCylinder.html |CircularCylinder|>, <CircularCylinder.html
% |CircularCylinder|>, <CircularShellCylinder.html
% |CircularShellCylinder|>, <EllipticCylinder.html |EllipticCylinder|>,
% <PolyognalCylinder.html |PolygonalCylinder|>, <Shape.html |Shape|>,
% <maxwell_run.html |maxwell_run|>

classdef SegmentalCylinder < GenericCylinder

	methods
        function this = SegmentalCylinder(normal_axis, height, center, radius, theta, d_theta, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');
			chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
				'"center" should be length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(radius, 'real') && radius > 0, '"radius" should be positive.');
			chkarg(istypesizeof(theta, 'real'), '"theta" should be real.');
			chkarg(istypesizeof(d_theta, 'real') && (d_theta < pi || d_theta > -pi), '"d_theta" should be real between -pi and pi.');
			
			sectoral = SectoralCylinder(normal_axis, height, center, radius, theta, d_theta, dl_max);
			
			thetas = [theta, theta + d_theta];
			
			mc_dir = mean([cos(thetas.') sin(thetas.')]);  % direction bisecting the segment
			mc_dir = mc_dir / norm(mc_dir);
			
			[h, v, n] = cycle(normal_axis);
			midpt = center([h, v]) + mean(radius .* [cos(thetas.') sin(thetas.')]);  % midpoint of chord

			% lsf2d() can handle rho = [p q] with column vectors p and q.  The
			% level set function is the one for a rectangle defined in the
			% (theta, radius) domain.
			function level = lsf2d(p, q)
				chkarg(istypeof(p, 'real'), '"p" should be array with real elements.');
				chkarg(istypeof(q, 'real'), '"q" should be array with real elements.');
				chkarg(isequal(size(p), size(q)), '"p" and "q" should have same size.');
				
				level1 = sectoral.lsf2d(p, q);  % level set function of sectoral cylinder

				lm = {p - midpt(Dir.h), q - midpt(Dir.v)};
				level2 = zeros(size(p));
				for d = Dir.elems
					level2 = level2 + lm{d} .* mc_dir(d);
				end
				
				level = min(level1, level2);  % intersection of regions defined by level1 and level2
			end
			
			lprim = cell(1, Axis.count);
			lprim{h} = radius * cos(thetas) + center(h);
			lprim{v} = radius * sin(thetas) + center(v);
			lprim{n} = [-height height]/2 + center(n);
			
			if sectoral.lsf_th(0) > 0
				lprim{h} = [lprim{h}, center(h) + radius];
			end
			if sectoral.lsf_th(pi/2) > 0
				lprim{v} = [lprim{v}, center(v) + radius];
			end
			if sectoral.lsf_th(pi) > 0
				lprim{h} = [lprim{h}, center(h) - radius];
			end
			if sectoral.lsf_th(3*pi/2) > 0
				lprim{v} = [lprim{v}, center(v) - radius];
			end
						
			if nargin < 7  % no dl_max
				super_args = {normal_axis, @lsf2d, lprim};
			else
				super_args = {normal_axis, @lsf2d, lprim, dl_max};
			end
			
			this = this@GenericCylinder(super_args{:});
		end
	end
end

