%% Sphere
% Concrete subclass of <Ellipsoid.html |Ellipsoid|> representing a sphere.

%%% Description
% |Sphere| represents the shape of a sphere.

%%% Construction
%  shape = Sphere(center, radius)
%  shape = Sphere(center, radius, dl_max)

% *Input Arguments*
%
% * |center|: center of the sphere in the format of |[x y z]|.
% * |radius|: radius of the sphere
% * |dl_max|: maximum grid size allowed in the sphere.  It can be either |[dx dy
% dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned, |dl_max
% = Inf| is used.

%%% Example
%   % Create an instance of Sphere.
%   shape = Sphere([0 0 0], 100);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <Ellipsoid.html |Ellipsoid|>, <EllipticCylinder.html |EllipticCylinder|>,
% <CircularCylinder.html |CircularCylinder|>, <Shape.html |Shape|>,
% <maxwell_run.html |maxwell_run|>

classdef Sphere < Ellipsoid
	% Sphere is a Shape for a sphere.

	methods
        function this = Sphere(center, radius, dl_max)
			chkarg(istypesizeof(radius, 'real') && radius > 0, '"radius" should be positive.');
			
			if nargin < 3  % no dl_max
				super_args = {center, [radius radius radius]};
			else
				super_args = {center, [radius radius radius], dl_max};
			end
			
			this = this@Ellipsoid(super_args{:});
		end
	end
end
