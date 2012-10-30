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
