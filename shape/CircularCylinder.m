classdef CircularCylinder < EllipticCylinder
	% CircularCylinder is a Shape for a cylinder whose cross section is a
	% circle.

	methods
        function this = CircularCylinder(normal_axis, center, radius, height, dl_max)
			chkarg(istypesizeof(radius, 'real') && radius > 0, '"radius" should be positive.');
			
			if nargin < 5
				super_args = {normal_axis, center, [radius radius], height};
			else
				super_args = {normal_axis, center, [radius radius], height, dl_max};
			end
			
			this = this@EllipticCylinder(super_args{:});
		end
	end
end

