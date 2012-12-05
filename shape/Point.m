%% Point
% Concrete subclass of <ZeroVolShape.html |Shape|> representing a point.

%%% Description
% |Point| represents the shape of a point.  It does not have a volume, but it is
% used to force a user-defined primary grid point in <maxwell_run.html
% |maxwell_run|>.

%%% Construction
%  shape = Point(location)

% *Input Arguments*
%
% * |location|: location of the point in the format of |[x y z]|.

%%% Example
%   % Create an instance of Point.
%   shape = Point(Axis.z, [0 0 0]);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <Plane.html |Plane|>, <Rectangle.html |Rectangle|>, <Line.html |Line|>,
% <ZeroVolShape.html |ZeroVolShape|>, <maxwell_run.html |maxwell_run|>

classdef Point < ZeroVolShape
	
	properties
		location  % [x, y, z]: location of point
	end
		
	methods
        function this = Point(location)
			chkarg(istypesizeof(location, 'real', [1, Axis.count]), ...
				'"intercept" should be length-%d row vector with real elements.', Axis.count);
			
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				
				N = size(r, 1);
				c = repmat(location, [N 1]);
				level = -max(abs(r - c), [], 2);
			end
			
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				lprim{w} = location(w);
			end
			
			super_args = {lprim, @lsf};  % dl_max is meaningless

			this = this@ZeroVolShape(super_args{:});
			this.location = location;
		end
	end	
end
