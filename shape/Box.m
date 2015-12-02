%% Box
% Concrete subclass of <Shape.html |Shape|> representing a 3D box.

%%% Description
% |Box| is the most basic shape.  Its sides are aligned with the axes of the
% Cartesian coordinate system.

%%% Construction
%  shape = Box(bound)
%  shape = Box(bound, dl_max)
% 
% *Input Arguments*
%
% * |bound|: bounds of the box in the format of |[xmin xmax; ymin ymax; zmin
% zmax]|.
% * |dl_max|: maximum grid size allowed in the box.  It can be either |[dx dy
% dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned, |dl_max
% = Inf| is used.

%%% Example
%   % Create an instance of Box.
%   shape = Box([-100 100; -50 50; 0 20]);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <Shape.html |Shape|>, <maxwell_run.html |maxwell_run|>

classdef Box < Shape
		
	methods
        function this = Box(bound, dl_max)
			% bound
			chkarg(istypesizeof(bound, 'real', [Axis.count, Sign.count]), ...
				'"bound" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
			s = diff(bound, 1, 2) ./ 2;  % semisides
			s = s.';  % row vector
			chkarg(all(s >= 0), '"bound" should have smaller lower bound than upper bound in all axes.');
			c = mean(bound, 2);  % center
			c = c.';  % row vector
			
			% The level set function is basically 1 - max(abs(r-c) ./ s, [], 2),
			% but it is vectorized, i.e., modified to handle r = [x y z] with
			% column vectors x, y, z.
			function level = lsf(x, y, z)
				chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
				chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
				chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
				chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');
				
				loc = {x, y, z};
				level = -Inf(size(x));
				
				for v = Axis.elems
					level = max(level, abs(loc{v}-c(v)) ./ s(v));
				end
				level = 1 - level;
			end
			
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				lprim{w} = bound(w,:);
			end
			
			if nargin < 2  % no dl_max
				super_args = {lprim, @lsf};
			else
				super_args = {lprim, @lsf, dl_max};
			end

			this = this@Shape(super_args{:});
		end
		
% 		function [n_dir, r_vol] = ndir_and_rvol(this, pixel)
% 			chkarg(istypesizeof(pixel, 'real', [Axis.count, Sign.count]), '"pixel" should be [xmin xmax; ymin ymax; zmin zmax].');
% 			center = mean(pixel, 2).';
% 			
% 			in = NaN(1, Axis.count);
% 			d = NaN(1, Axis.count);
% 			for w = Axis.elems
% 				[in(w), d(w)] = this.interval(w).contains(center(w));
% 			end
% 			chkarg((in(Axis.x) && in(Axis.y)) || (in(Axis.y) && in(Axis.z)) || (in(Axis.z) && in(Axis.x)), 'center of "pixel" should be contained the box at least in two directions.');
% 
% 			if all(in)
% 				[~, i] = min(d);
% 			else  % only one element of "in" is false
% 				i = find(~in);
% 			end
% 			
% 			w = Axis.elems(i);
% 			n_dir = zeros(1,Axis.count);
% 			n_dir(w) = 1;
% 			dc = d(w);
% 			L = diff(pixel(w,:));
% 			chkarg(L, '"pixel" should be [xmin xmax; ymin ymax; zmin zmax], where wmin < wmax.');
% 			if dc >= L/2
% 				r_vol = 0;
% 			elseif in(w)
% 				r_vol = (L/2 + dc) / L;
% 			else  % in(w) == false
% 				r_vol = (L/2 - dc) / L;
% 			end
% 		end
	end	
end
