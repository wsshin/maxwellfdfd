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
			
			function level = lsf(x, y, z)
				chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
				chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
				chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
				chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');
				
				loc = {x, y, z};
				level = -Inf(size(x));
				for v = Axis.elems
					level = max(level, abs(loc{v} - location(v)));
				end
				level = -level;
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
