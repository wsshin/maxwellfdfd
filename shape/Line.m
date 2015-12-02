%% Line
% Concrete subclass of <ZeroVolShape.html |Shape|> representing a plane.

%%% Description
% |Line| represents the shape of a line.  It does not have a volume, but it is
% used to force a user-defined primary grid line in <maxwell_run.html
% |maxwell_run|>.

%%% Construction
%  shape = Line(axis, intercept)
%  shape = Line(axis, intercept, dl_max)

% *Input Arguments*
%
% * |axis|: axis of the line.  It should be one of |Axis.x|, |Axis.y|, |Axis.z|.
% * |intercept|: location of the line in the plane normal to the line.  For
% |axis == Axis.y|, it is in the format of |[z x]|.
% * |dl_max|: maximum grid size allowed in the plane.  It can be either |[dx dy
% dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned, |dl_max
% = Inf| is used.  In the directions normal to the |axis| direction, |dl_max| is
% meaningless.

%%% Example
%   % Create an instance of Line.
%   shape = Line(Axis.z, [0 100]);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <Plane.html |Plane|>, <Rectangle.html |Rectangle|>, <Point.html |Point|>,
% <ZeroVolShape.html |ZeroVolShape|>, <maxwell_run.html |maxwell_run|>

classdef Line < ZeroVolShape
	
	properties
		axis  % direction of line
		intercept  % [h, v]: location of intercept in transverse plane
	end
		
	methods
        function this = Line(axis, intercept, dl_max)
			chkarg(istypesizeof(axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real', [1, Dir.count]), ...
				'"intercept" should be length-%d row vector with real elements.', Dir.count);
			
			[h, v] = cycle(axis);  % h-axis, v-axis
			
			function level = lsf(x, y, z)
				chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
				chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
				chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
				chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');
				
				loc = {x, y, z};
				level = -max(abs(loc{h} - intercept(Dir.h)), abs(loc{v} - intercept(Dir.v)));
			end
			
			lprim = cell(1, Axis.count);
			lprim{axis} = [-Inf Inf];
			lprim{h} = intercept(Dir.h);
			lprim{v} = intercept(Dir.v);
			
			if nargin < 3  % no dl_max
				super_args = {lprim, @lsf};
			else
				dl_max = expand2row(dl_max, Axis.count);
				dl_max(h) = Inf;  dl_max(v) = Inf;  % dl_max is meaningful only in direction of line
				super_args = {lprim, @lsf, dl_max};
			end

			this = this@ZeroVolShape(super_args{:});
			this.axis = axis;
			this.intercept = intercept;
		end
	end	
end
