%% ZeroVolShape
% Subclass of <Shape.html |Shape|> representing a shape without volume.

%%% Description
% |ZeroVolShape| is the superclass of shapes without volume such as a plane,
% point, line.  |ZeroVolShape| is designed to help visualization of shapes
% without volume, because such shapes are hard to visualize.

%%% See Also
% <Plane.html |Plane|>, <Rectangle.html |Rectangle|>, <Line.html |Line|>,
% <Point.html |Point|>, <maxwell_run.html |maxwell_run|>

classdef ZeroVolShape < Shape

% 	methods (Abstract)
% 		draw2d(this, axes_handle, normal_axis, intercept)
% 		draw3d(this, axes_handle)
% 	end
	
	methods
		function this = ZeroVolShape(lprim_cell, lsf, dl_max)
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
				super_args = {lprim_cell, @lsf_zv};
			else
				super_args = {lprim_cell, @lsf_zv, dl_max};
			end
			this = this@Shape(super_args{:});
		end
	end
end
