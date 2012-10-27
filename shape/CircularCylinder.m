classdef CircularCylinder < EllipticCylinder
	% CircularCylinder is a Shape for a cylinder whose cross section is a
	% circle.

	methods
        function this = CircularCylinder(normal_axis, center, radius, height, dl_max, dl_boundary)
			chkarg(istypesizeof(radius, 'real') && radius > 0, '"radius" should be positive.');
			
			if nargin < 5
				super_args = {normal_axis, center, [radius radius], height};
			elseif nargin < 6
				super_args = {normal_axis, center, [radius radius], height, dl_max};
			else
				super_args = {normal_axis, center, [radius radius], height, dl_max, dl_boundary};
			end
			
			this = this@EllipticCylinder(super_args{:});
		end
	end
end

