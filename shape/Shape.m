classdef Shape < handle & matlab.mixin.Heterogeneous
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
	
	properties (Access = public)
		n_subpxls  % number of uniform sampling points in each Cartesian direction for subpixel smoothing
	end
	
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
			
			this.n_subpxls = 10;
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
		
		function truth = circumbox_contains(this, x, y, z)
			chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
			chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
			chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
			chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');

			truth = true(size(x));
			loc = {x, y, z};
			for w = Axis.elems
				truth = truth & this.comp(w).contains(loc{w});  % &: elementwise AND operator
			end
		end
		
		function truth = contains(this, x, y, z)
			chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
			chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
			chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
			chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');

			truth = this.lsf(x, y, z) >= 0;
			
% 			% Check if circumbox is correctly set.
% 			if truth
% 				chkarg(this.circumbox_contains(point), ...
% 					'circumbox of this shape does not contain the shape.');
% 			end
		end
		
		function [rvol, ndir] = smoothing_params(this, box)
			chkarg(istypesizeof(box, 'real', [Axis.count, Sign.count]), ...
				'"box" should be %d-by-%d array with real elements.', Axis.count, Sign.count);

			l_probe = cell(1, Axis.count);	
			dl_sub = diff(box.') ./ this.n_subpxls;
			N_subpxls = this.n_subpxls^3;
			
			for w = Axis.elems
				bnd = box(w, :);  % interval in w-direction
				
				% Below, probing points are at the centers of subpixels.  One
				% extra subpixel beyond each boundary of the box is considered
				% to discard single-sided differences in gradient().
				l_probe{w} = linspace(bnd(Sign.n) - dl_sub(w)/2, bnd(Sign.p) + dl_sub(w)/2, this.n_subpxls + 2);
			end
			
			[X, Y, Z] = ndgrid(l_probe{:});
			F = this.lsf(X, Y, Z);
			
			Fint = F(2:end-1, 2:end-1, 2:end-1);  % values at interior points of box
			is_contained = (Fint > 0);
			at_interface = (Fint == 0);
			rvol = sum(double(is_contained(:))) + 0.5 * sum(double(at_interface(:)));
			rvol = rvol / N_subpxls;
			
			[Fy, Fx, Fz] = gradient(F, dl_sub(Axis.y), dl_sub(Axis.x), dl_sub(Axis.z));  % x and y are swapped due to MATLAB's definiton of gradient
			Fx = Fx(2:end-1, 2:end-1, 2:end-1);
			Fy = Fy(2:end-1, 2:end-1, 2:end-1);
			Fz = Fz(2:end-1, 2:end-1, 2:end-1);
			
			gradF = [Fx(:), Fy(:), Fz(:)];
			normF = sqrt(sum(gradF.^2, 2));
			gradF = bsxfun(@rdivide, gradF, normF);
			
			ndir = sum(gradF) ./ N_subpxls;
			ndir = -ndir ./ norm(ndir);
		end
	end
end
