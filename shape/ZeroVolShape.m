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
			function level = lsf_zv(x, y, z, force_draw)
				if nargin < 4  % no force_draw
					force_draw = false;
				end
				
				level = lsf(x, y, z);
				level_copy = level(:);
				if force_draw && all(level_copy <= 0)
					% If the level set function is negative, nothing is drawn,
					% so shift the function upward.
					max_level = max(level_copy);
					max_level = max(level_copy(level_copy<max_level));  % second largest
% 					max_level = max(level_copy(level_copy<max_level));  % third largest
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
