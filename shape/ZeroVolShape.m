classdef ZeroVolShape < Shape
	% ZeroVolShape is a superclass for all shapes with zero volume.

% 	methods (Abstract)
% 		draw2d(this, axes_handle, normal_axis, intercept)
% 		draw3d(this, axes_handle)
% 	end
	
	methods
		function this = ZeroVolShape(circumbox, lsf, dl_max)
			function level = lsf_zv(r, force_draw)
				if nargin < 2  % no force_draw
					force_draw = false;
				end
				
				level = lsf(r);
				if force_draw && all(level <= 0)
					% Make sure to "level" reaches zero; otherwise the plane is not
					% drawn.
					max_level = max(level);
					level = level - max_level;
				end
			end

			if nargin < 3  % no dl_max
				super_args = {circumbox, @lsf_zv};
			else
				super_args = {circumbox, @lsf_zv, dl_max};
			end
			this = this@Shape(super_args{:});
		end
	end
end
