classdef Plane < ZeroVolShape
	% Plane is a Shape for a plane.
	
	properties
		normal_axis
		intercept
	end
		
	methods
        function this = Plane(normal_axis, intercept, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			
			bound = [-Inf Inf; -Inf Inf; -Inf Inf];
			bound(normal_axis, :) = [intercept intercept];
			
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				r_normal = r(:, normal_axis);
				level = -abs(r_normal - intercept);				
			end
			
			if nargin < 3  % no dl_max
				super_args = {bound, @lsf};
			else
				dl_max = expand2row(dl_max, Axis.count);
				dl_max(normal_axis) = Inf;  % dl_max is meaningless in normal direction
				super_args = {bound, @lsf, dl_max};
			end

			this = this@ZeroVolShape(super_args{:});
			this.normal_axis = normal_axis;
			this.intercept = intercept;
		end
	end	
end
