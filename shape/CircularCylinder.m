%% CircularCylinder
% Concrete subclass of <EllipticCylinder.html |EllipticCylinder|> representing a
% cylinder with a circular cross section.

%%% Description
% |CircularCylinder| represents the shape of a circular cylinder.  The axis of
% the cylinder should be aligned with one of the axes of the Cartesian
% coordinate system.

%%% Construction
%  shape = CircularCylinder(normal_axis, height, center, radius)
%  shape = CircularCylinder(normal_axis, height, center, radius, dl_max)
% 
% *Input Arguments*
%
% * |normal_axis|: axis of the cylinder.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |height|: size of the cylinder along its axis.
% * |center|: center of the cylinder in the format of |[x y z]|.
% * |radius|: radius of the circle
% * |dl_max|: maximum grid size allowed in the cylinder.  It can be either |[dx
% dy dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned,
% |dl_max = Inf| is used.

%%% Example
%   % Create an instance of CircularCylinder.
%   shape = CircularCylinder(Axis.z, 100, [0 0 50], 50);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <CircularShellCylinder.html |CircularShellCylinder|>, <EllipticCylinder.html
% |EllipticCylinder|>, <SectoralCylinder.html |SectoralCylinder|>,
% <PolyognalCylinder.html |PolygonalCylinder|>, <Shape.html |Shape|>,
% <maxwell_run.html |maxwell_run|>

classdef CircularCylinder < EllipticCylinder
	% CircularCylinder is a Shape for a cylinder whose cross section is a
	% circle.

	methods
        function this = CircularCylinder(normal_axis, height, center, radius, dl_max)
			chkarg(istypesizeof(radius, 'real') && radius > 0, '"radius" should be positive.');
			
			if nargin < 5
				super_args = {normal_axis, height, center, [radius radius]};
			else
				super_args = {normal_axis, height, center, [radius radius], dl_max};
			end
			
			this = this@EllipticCylinder(super_args{:});
		end
	end
end

