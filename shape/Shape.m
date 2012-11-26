classdef Shape < matlab.mixin.Heterogeneous
	% Shape is a superclass for all shapes.
	
	properties (SetAccess = immutable)
		interval  % [interval_x, interval_y, interval_z]
		
		% Below, "lsf" is a level set function.  "lsf" takes a position r =
		% [x,y,z] as an argument and returns a real number.  For [x,y,z] that is
		% in the interior, boundary, and exterior of the shape, lsf(r) >0, ==0,
		% and <0. See http://en.wikipedia.org/wiki/Level_set_method for more
		% details. In addition, "lsf" must be vectorized, i.e., handle column
		% vectors x, y, z.
		lsf  % level set function
	end
	
	properties (Dependent, SetAccess = immutable)
		lprim  % {xprim_array, yprim_array, zprim_array}: primary grid plane that this shape defines
		bound  % [xmin xmax; ymin ymax; zmin zmax]: circumbox of this shape
		cb_center  % center of circumbox
		L  % [Lx, Ly, Lz]: range of circumbox
		dl_max  % [dx_max, dy_max, dz_max]: maximum grid cell size in this shape
	end		
	
% 	methods (Abstract)  % no abstract methods allowed for periodize_object()
% 		[n_dir, r_vol] = ndir_and_rvol(this, pixel)  % for subpixel smoothing
% 	end
	
	methods
		function this = Shape(lprim_cell, lsf, dl_max)
			% lprim_cell
			chkarg(istypesizeof(lprim_cell, 'realcell', [1 Axis.count], [1 0]), ...
				'"lprim_cell" is length-%d row cell array whose each element is row vector with real elements.', Axis.count);
			
			% level set function
			chkarg(istypeof(lsf, 'function_handle'), '"lsf" should be function handle.');
			this.lsf = lsf;
			
			% dl_max
			if nargin < 3  % no dl_max
				dl_max = Inf;
			end
			chkarg(istypeof(dl_max, 'real') && all(dl_max > 0), 'element of "dl_max" should be positive.');
			chkarg(isexpandable2row(dl_max, Axis.count), ...
				'"dl_max" should be scalar or length-%d vector.', Axis.count);
			dl_max = expand2row(dl_max, Axis.count);

			this.interval = Interval.empty();
			for w = Axis.elems
				this.interval(w) = Interval(lprim_cell{w}, dl_max(w));
			end
		end
		
		function lprim = get.lprim(this)
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				lprim{w} = this.interval(w).lprim;
			end
		end
		
		function bound = get.bound(this)
			bound = NaN(Axis.count, Sign.count);
			for w = Axis.elems
				bound(w,:) = this.interval(w).bound;
			end
		end
		
		function center = get.cb_center(this)
			center = mean(this.bound, 2);
			center = center.';
		end
				
		function L = get.L(this)
			L = NaN(1, Axis.count);
			for w = Axis.elems
				L(w) = this.interval(w).L;
			end
		end
		
		function dl_max = get.dl_max(this)
			dl_max = NaN(1, Axis.count);
			for w = Axis.elems
				dl_max(w) = this.interval(w).dl_max;
			end
		end
		
		function truth = circumbox_contains(this, point, axes)
			chkarg(istypesizeof(point, 'real', [0 Axis.count]), ...
				'"point" should be matrix with %d columns with real elements.', Axis.count);
			
			if nargin < 3  % no axes
				axes = Axis.elems;
			end
			chkarg(istypesizeof(axes, 'Axis', [1 0]) && length(axes) <= Axis.count, ...
				'"axes" should be length-%d or shorter row vector of instances of Axis.', Axis.count);

			n = size(point, 1);
			truth = true(n, 1);
			for w = axes
				truth = truth & this.interval(w).contains(point(:,w));
			end
		end
		
		function truth = contains(this, point)
			chkarg(istypesizeof(point, 'real', [0 Axis.count]), ...
				'"point" should be matrix with %d columns with real elements.', Axis.count);
			truth = this.lsf(point) >= 0;
			
% 			% Check if circumbox is correctly set.
% 			if truth
% 				chkarg(this.circumbox_contains(point), ...
% 					'circumbox of this shape does not contain the shape.');
% 			end
		end
	end
end
