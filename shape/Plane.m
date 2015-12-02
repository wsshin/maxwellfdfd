%% Plane
% Concrete subclass of <ZeroVolShape.html |Shape|> representing a plane.

%%% Description
% |Plane| represents the shape of a plane.  It does not have a volume, but it is
% used to force a user-defined primary grid plane in <maxwell_run.html
% |maxwell_run|>.

%%% Construction
%  shape = Plane(normal_axis, intercept)
%  shape = Plane(normal_axis, intercept, dl_max)

% *Input Arguments*
%
% * |normal_axis|: axis normal to the plane.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |intercept|: location of the plane in the |normal_axis| direction.
% * |dl_max|: maximum grid size allowed in the plane.  It can be either |[dx dy
% dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned, |dl_max
% = Inf| is used.  In the |normal_axis| direction, |dl_max| is meaningless.

%%% Example
%   % Create an instance of Plane.
%   shape = Plane(Axis.y, 100);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <Rectangle.html |Rectangle|>, <Line.html |Line|>, <Point.html |Point|>,
% <ZeroVolShape.html |ZeroVolShape|>, <maxwell_run.html |maxwell_run|>

classdef Plane < ZeroVolShape
	
	properties
		normal_axis
		intercept
	end
		
	methods
        function this = Plane(normal_axis, intercept, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			
			function level = lsf(x, y, z)
				chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
				chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
				chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
				chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');
				
				loc = {x, y, z};
				level = -abs(loc{normal_axis} - intercept);				
			end
			
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				if w == normal_axis
					lprim{w} = intercept;
				else
					lprim{w} = [-Inf Inf];
				end
			end
			
			if nargin < 3  % no dl_max
				super_args = {lprim, @lsf};
			else
				dl_max = expand2row(dl_max, Axis.count);
				dl_max(normal_axis) = Inf;  % dl_max is meaningless in normal direction
				super_args = {lprim, @lsf, dl_max};
			end

			this = this@ZeroVolShape(super_args{:});
			this.normal_axis = normal_axis;
			this.intercept = intercept;
		end
	end	
end
