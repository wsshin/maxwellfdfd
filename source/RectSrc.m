%% RectSrc
% Concrete subclass of <Source.html |Source|> representing a constant electric
% dipole distribution over a rectangular patch.

%%% Description
% |RectSrc| is similar to <PlaneSrc.html |PlaneSrc|>, but it places constant
% electric dipoles on a rectangular patch rather than an entire plane.

%%% Construction
%  src = RectSrc(normal_axis, intercept, rect, polarization)
%  src = RectSrc(normal_axis, intercept, rect, polarization, K)
% 
% *Input Arguments*
%
% * |normal_axis|: direction normal to the rectangular patch.  It should be one
% of |Axis.x|, |Axis.y|, |Axis.z|.
% * |intercept|: location of the rectangular patch in the |normal_axis|
% direction.
% * |rect|: description of the rectangular patch in the format of |[hmin hmax;
% vmin vmax]|.  If |normal_axiz == Axis.y|, then it is |[zmin zmax; xmin xmax]|.
% * |polarization|: direction of the dipoles distributed on the rectangular
% patch.  It should be one of |Axis.x|, |Axis.y|, |Axis.z| that is orthogonal to
% |normal_axis|.
% * |K|: amplitude of the surface current density that the dipoles drive.

%%% Note
% In the finite-difference grid, |RectSrc| excites dipoles at the _E_-field
% points.  This poses a condition on |intercept| argument in the constructor:
% |intercept| should be at a dual grid point in the |normal_axis| direction.
% Therefore, make sure that |intercept| does not overlap with the locations of
% the vertices of <Shape.html |Shape|> in the |normal_axis| direction; otherwise
% dynamic grid generation in <moxwell_run.html |maxwell_run|> will fail.

%%% Example
%   % Create an instance of PointSrc.
%   src =  RectSrc(Axis.y, -2000, [0 10; 20 180], Axis.x);  % y = -2000 should not be primary grid point
%
%   % Use the constructed src in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'SRCJ', src);

%%% See Also
% <PlaneSrc.html |PlaneSrc|>, <maxwell_run.html |maxwell_run|>

classdef RectSrc < Source
	
	properties (SetAccess = immutable)
		normal_axis  % plane normal axis: one of Axis.x, Axis.y, Axis.z
		intercept  % intercept between plane and normal axis
		rect  % rectangular patch
		polarization  % one of Axis.x, Axis.y, Axis.z
		K  % sheet current density
	end
	
	methods
		function this = RectSrc(normal_axis, intercept, rect, polarization, K)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			chkarg(istypesizeof(rect, 'real', [Dir.count Sign.count]), ...
				'"rect" should be %d-by-%d array with real elements.', Dir.count, Sign.count);
			chkarg(istypesizeof(polarization, 'Axis') && polarization ~= normal_axis , ...
				'"polarization" should be instance of Axis and different from "normal_axis".');
			
			if nargin < 5  % no K
				K = 1.0;
			end
			chkarg(istypesizeof(K, 'complex'), '"K" should be complex.');

			lgrid = cell(1, Axis.count);
			laltgrid = cell(1, Axis.count);
			lgrid{normal_axis} = intercept;
			
			rectangle = Rectangle(normal_axis, intercept, rect);
			this = this@Source(lgrid, laltgrid, rectangle);
			
			this.normal_axis = normal_axis;
			this.intercept = intercept;
			this.rect = rect;
			this.polarization = polarization;
			this.K = K;
		end
		
		function [index_cell, JMw_patch] = generate_kernel(this, w_axis, grid3d)
			index_cell = cell(1, Axis.count);
			if w_axis ~= this.polarization
				JMw_patch = [];
			else  % w_axis == this.polarization
				n = this.normal_axis;
				p = w_axis;  % polarization axis
% 				q = setdiff(Axis.elems, [n, p]);
				
				g = this.gt;
				ind_n = ind_for_loc(this.intercept, n, g, grid3d);
				dn = grid3d.dl{n,g}(ind_n);
				JM = this.K / dn;  % for t normal to both n and K, K*dt = (current through dn*dt) = J * (dn*dt)

				index_cell{n} = ind_n;
				r_overlap = cell(1, Dir.count);  % ratios of intervals overlapping with this rectangle
				grid2d = Grid2d(grid3d, n);
				for d = Dir.elems
					if grid2d.axis(d) == p  % overlap in polarization direction
						g = alter(this.gt);
					else  % overlap in in-plane direction normal to polarization
						g = this.gt;
					end
					[index_cell{grid2d.axis(d)}, r_overlap{d}] = this.get_ind_r(d, g, grid2d);
				end
				
				JMw_patch = JM .* (r_overlap{Dir.h}.' * r_overlap{Dir.v});
				JMw_patch = ipermute(JMw_patch, int([grid2d.axis(Dir.h) grid2d.axis(Dir.v) n]));
			end
		end		
	end
	
	methods (Access = private)
		function [ind, r_overlap] = get_ind_r(this, dir, g, grid2d)
			chkarg(istypesizeof(dir, 'Dir'), '"dir" should be instance of Dir.');
			chkarg(istypesizeof(g, 'GT'), '"g" should be instance of GT.');  % grid kind on which J will be assigned in dir-direction
			chkarg(istypesizeof(grid2d, 'Grid2d'), '"grid2d" should be instance of Grid2d.');
			
			interval = this.rect(dir, :);
			grid1d = grid2d.comp(dir);
			
			exception = MException.empty(0,1);
			if diff(interval) > diff(grid1d.bound)
				exception = MException('Maxwell:srcAssign', 'this %s should be narrower than grid in %s-axis.', ...
								class(this), char(grid2d.axis(dir)));
			end
			if ~any(grid1d.contains(interval.'))
				exception = MException('Maxwell:srcAssign', 'this %s should have overlap with grid in %s-axis.', ...
								class(this), char(grid2d.axis(dir)));
			end
			
			if ~isempty(exception)
				throw(exception);
			end

			% Adjust interval as if the boundary condition is periodic;
			% nonperiodic BC will be handled later.
			interval_adjusted = false;
			if ~grid1d.contains(interval(Sign.n))
				interval(Sign.n) = interval(Sign.n) + diff(grid1d.bound);
				interval_adjusted = true;
				drop_n = true;  % if BC is nonperiodic, drop negative portion
			end
			
			if ~grid1d.contains(interval(Sign.p)) % || (interval(Sign.p) == grid1d.bound(Sign.p) && grid1d.bc(Sign.n) == BC.p)
				interval(Sign.p) = interval(Sign.p) - diff(grid1d.bound);
				interval_adjusted = true;
				drop_n = false;  % if BC is nonperiodic, drop positive portion
			end
			
			% Suppose that g == GT.prim.  Then, the purpose below is to find
			% the primary grid indices for which J are assigned.  Note that J at
			% each primary grid extends to the dual grids around it.
			ind_n = find(grid1d.lall{alter(g)} <= interval(Sign.n) , 1, 'last');  % grid1d.lg rather than grid1d.l
			ind_p = find(grid1d.lall{alter(g)} >= interval(Sign.p) , 1, 'first') - 1;  % grid1d.lg rather than grid1d.l
			if isempty(ind_p)
				exception = MException('Maxwell:srcAssign', 'this %s extends beyond last %s grid point in %s-axis.', ...
								class(this), char(alter(g)), char(grid2d.axis(dir)));
				throw(exception);
			end
			
			if ind_p == grid1d.N+1
				ind_p = 1;
				interval_adjusted = true;
			end
			
			if ind_n < ind_p
				ind = ind_n:ind_p;
			elseif ind_n == ind_p  % can happen by drop_n assignment or mod(ind_p - 1, ...)
				if interval_adjusted % && grid1d.bc(Sign.n) == BC.p
					ind = 1:grid1d.N;
				else
					ind = ind_n;
				end
			else  % ind_n > ind_p
				ind = [1:ind_p, ind_n:grid1d.N];
			end
			
			r_overlap = ones(1, length(ind));
			if ind_n < ind_p  % ind = ind_n:ind_p;
				r_overlap(1) = (grid1d.lg{alter(g)}(ind_n+1) - interval(Sign.n)) / grid1d.dl{g}(ind_n);
				r_overlap(end) = (interval(Sign.p) - grid1d.lg{alter(g)}(ind_p)) / grid1d.dl{g}(ind_p);
			elseif ind_n == ind_p
				if interval_adjusted  % ind = 1:grid1d.N;
					r_overlap(ind_n) = ...
						(grid1d.lg{alter(g)}(ind_n+1) - interval(Sign.n) ...
						+ interval(Sign.p) - grid1d.lg{alter(g)}(ind_p)) ...
						/ grid1d.dl{g}(ind_n);
				else  % ind = ind_n;
					r_overlap(1) = diff(interval) / grid1d.dl{g}(ind_n);
				end
			else  % ind_n > ind_p: ind = [1:ind_p, ind_n:grid1d.N];
				r_overlap(ind_p) = (interval(Sign.p) - grid1d.lg{alter(g)}(ind_p)) / grid1d.dl{g}(ind_p);
				r_overlap(ind_p+1) = (grid1d.lg{alter(g)}(ind_n+1) - interval(Sign.n)) / grid1d.dl{g}(ind_n);
			end
			
			if ind_n > ind_p && grid1d.bc(Sign.n) ~= BC.p  % cannot handle case with interval_adjusted == true and ind_n == ind_p
				if drop_n
					ind = ind(1:ind_p);
					r_overlap = r_overlap(1:ind_p);
				else
					ind = ind(ind_p+1:end);
					r_overlap = r_overlap(ind_p+1:end);
				end
			end
		end
	end
end

