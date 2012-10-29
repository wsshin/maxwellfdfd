classdef Plane < Shape
	% Plane is a Shape for a plane.
	
	properties
		normal_axis
		intercept
	end
		
	methods
        function this = Plane(normal_axis, intercept, dl_max, dl_boundary)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			
			bound = [-Inf Inf; -Inf Inf; -Inf Inf];
			bound(normal_axis, :) = [intercept intercept];
			
			function level = lsf(r, force_draw)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				if nargin < 2  % no force_draw
					force_draw = false;
				end
				r_normal = r(:, normal_axis);
				level = -abs(r_normal - intercept);
				
				if force_draw
					% Make sure to "level" reaches zero; otherwise the plane is not
					% drawn.
					max_level = max(level);
					level = level - max_level;
				end
			end
			
			if nargin < 3  % no dl_max
				super_args = {bound, @lsf};
			elseif nargin < 4  % no dl_boundary
				super_args = {bound, @lsf, dl_max};
			else
				super_args = {bound, @lsf, dl_max, dl_boundary};
			end

			this = this@Shape(super_args{:});
			this.normal_axis = normal_axis;
			this.intercept = intercept;
		end
	end	
end
