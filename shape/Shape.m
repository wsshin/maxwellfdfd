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
		lsf_orig  % level set function
	end
	
	properties (Dependent, SetAccess = immutable)
		lprim  % {xprim_array, yprim_array, zprim_array}: primary grid plane that this shape defines
		bound  % [xmin xmax; ymin ymax; zmin zmax]: circumbox of this shape
		cb_center  % center of circumbox
		L  % [Lx, Ly, Lz]: range of circumbox
		dl_max  % [dx_max, dy_max, dz_max]: maximum grid cell size in this shape
		lsf  % intersection of lsf_orig and domain.lsf
	end
	
	properties (Access = public)
		n_subpxls  % number of uniform sampling points in each Cartesian direction for subpixel smoothing
	end
	
	properties (Access = private)
		domain_  % instance of Domain
	end
	
	properties (Dependent)
		domain
	end
	
% 	methods (Sealed)
% 		function truth = ne(this, another)
% 			chkarg(isequal(size(this), size(another)) || numel(another) == 1 , '"this" and "another" should have same size, or another is scalar.');
% 			dims = size(this);
% 			n = numel(this);
% 			this = this(:);
% 			another = another(:);
% 			truth = true(n,1);
% 			if numel(another) == 1
% 				for i = 1:n
% 					truth(i)= ~(this(i) == another);
% 				end
% 			else
% 				for i = 1:n
% 					truth(i)= ~(this(i) == another(i));
% 				end
% 			end
% 			truth = reshape(truth, dims);
% 		end
% 	end
	
	methods
		function this = Shape(lprim_cell, lsf, dl_max)
			chkarg(istypesizeof(lprim_cell, 'realcell', [1 Axis.count], [1 0]), ...
				'"lprim_cell" is length-%d row cell array whose each element is row vector with real elements.', Axis.count);
			chkarg(istypeof(lsf, 'function_handle'), '"lsf" should be function handle.');
			
			this.lsf_orig = lsf;
			
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
			this.domain_ = Domain.empty();
		end
		
		function domain = get.domain(this)
			domain = this.domain_;
		end
		
		function set.domain(this, domain)
			chkarg(istypesizeof(domain, 'Domain'), '"domain" should be instance of Domain.');
			if ~isempty(this.domain)
				exception = MException('Maxwell:stateChk', 'Wrong state: "domain" is already set.');
				throwAsCaller(exception);
			end
			
			this.domain_ = domain;
		end
		
		function lsf = get.lsf(this)
			if isempty(this.domain)
				lsf = this.lsf_orig;
			else
				lsf = @(x, y, z) min(this.lsf_orig(x, y, z), this.domain.lsf_orig(x, y, z));
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
			
			if ~isempty(this.domain)
				truth = truth & this.domain.contains(x, y, z);
			end
		end
		
		function truth = contains(this, x, y, z)
			chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
			chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
			chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
			chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');

			% Below, the inequality must be >= 0 rather than >0.  Consider Ey at
			% the -x-boundary near the top-left corner of a box.
			truth = this.lsf(x, y, z) >= 0;
			
% 			% Check if circumbox is correctly set.
% 			if truth
% 				chkarg(this.circumbox_contains(point), ...
% 					'circumbox of this shape does not contain the shape.');
% 			end
		end
		
		function shape = flip_domain(this, axis, sign)
			if isempty(this.domain)
				exception = MException('Maxwell:stateChk', 'Wrong state: "domain" is empty.');
				throwAsCaller(exception);
			end
			
			chkarg(istypeof(axis, 'Axis'), '"axis" should be instance of Axis.');
			chkarg(istypeof(sign, 'Sign'), '"sign" should be instance of Sign.');

			w = axis;
			db_flipped = this.domain.bound;  % domain bound flipped
			w_flip = db_flipped(w, sign);
			w_db = db_flipped(w, :);
			db_flipped(w, :) = 2*w_flip - fliplp(w_db);
			domain_flipped = Domain(db_flipped, this.domain.dl_max);
			
			lprim_flipped = this.lprim;
			lprim_flipped{w} = 2*w_flip - fliplp(lprim_flipped{w});
			
			function level = lsf_flipped(x, y, z)
				loc_flipped = {x, y, z};
				loc_flipped{w} = 2*w_flip - loc_flipped{w};
				
				level = this.lsf_orig(loc_flipped{:});
			end
			
			shape = Shape(lprim_flipped, lsf_flipped, this.dl_max);
			shape.domain = domain_flipped;
		end
		
		function shape = shift_domain(this, axis, sign)
			if isempty(this.domain)
				exception = MException('Maxwell:stateChk', 'Wrong state: "domain" is empty.');
				throwAsCaller(exception);
			end
			
			chkarg(istypeof(axis, 'Axis'), '"axis" should be instance of Axis.');
			chkarg(istypeof(sign, 'Sign'), '"sign" should be instance of Sign.');

			w = axis;
			if sign == Sign.n
				s = -1.0;
			else  % sign == Sign.p
				s = 1.0;
			end

			db_shifted = this.domain.bound;  % domain bound shifted
			w_db = db_shifted(w, :);
			w_shift = s * diff(w_db);
			db_shifted(w, :) = w_db + w_shift;
			domain_shifted = Domain(db_shifted, this.domain.dl_max);

			
			lprim_shifted = this.lprim;
			lprim_shifted{w} = lprim_shifted{w} + w_shift;

			function level = lsf_shifted(x, y, z)
				loc_shifted = {x, y, z};
				loc_shifted{w} = loc_shifted{w} - w_shift;
				
				level = this.lsf_orig(loc_shifted{:});
			end

			shape = Shape(lprim_shifted, @lsf_shifted, this.dl_max);
			shape.domain = domain_shifted;
		end
		
		function rvol = fill_factor(this, voxel)
			chkarg(istypesizeof(voxel, 'real', [Axis.count, Sign.count]), ...
				'"voxel" should be %d-by-%d array with real elements.', Axis.count, Sign.count);

			l_probe = cell(1, Axis.count);	
			dl_sub = diff(voxel.') ./ this.n_subpxls;
			N_subpxls = this.n_subpxls^3;
			
			for w = Axis.elems
				bnd = voxel(w, :);  % interval in w-direction
				l_probe{w} = linspace(bnd(Sign.n) + dl_sub(w)/2, bnd(Sign.p) - dl_sub(w)/2, this.n_subpxls);
			end
			
			[X, Y, Z] = ndgrid(l_probe{:});
			F = this.lsf(X, Y, Z);
			
			is_contained = (F > 0);
			at_interface = (F == 0);
			rvol = sum(double(is_contained(:))) + 0.5 * sum(double(at_interface(:)));
			rvol = rvol / N_subpxls;
		end
		
		function ndir = outward_normal(this, voxel)
			chkarg(istypesizeof(voxel, 'real', [Axis.count, Sign.count]), ...
				'"voxel" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
			
			ndir = this.outward_normal_center(voxel);
			if norm(ndir) == 0  % possile if voxel center is at crease of lsf in union shape
				ndir = this.outward_normal_edge_avg(voxel);
			end
		end
		
		function ndir = outward_normal_edge_avg(this, voxel)
			chkarg(istypesizeof(voxel, 'real', [Axis.count, Sign.count]), ...
				'"voxel" should be %d-by-%d array with real elements.', Axis.count, Sign.count);

			dl = diff(voxel.');
			l = mat2cell(voxel.', Sign.count, ones(1, Axis.count));
			[X, Y, Z] = ndgrid(l{:});
			F = this.lsf(X, Y, Z);
			
			ndir = NaN(1, Axis.count);
			for w = Axis.elems
				dFw = -diff(F, 1, int(w));
				ndir(w) = mean(dFw(:)) / dl(w);
			end
			
			norm_n = norm(ndir);
			if norm_n < 10*eps
				ndir = [0 0 0];
			else
				ndir = ndir ./ norm_n;
			end
		end
		
		function ndir = outward_normal_center(this, box)
			chkarg(istypesizeof(box, 'real', [Axis.count, Sign.count]), ...
				'"box" should be %d-by-%d array with real elements.', Axis.count, Sign.count);

			c = mean(box.');
			dl = diff(box.');
			l = mat2cell(repmat(c, Sign.count*Axis.count, 1), Sign.count*Axis.count, ones(1, Axis.count));
			for w = Axis.elems
				for s = Sign.elems
					l{w}(Sign.count * subsindex(w) + int(s)) = box(w, s);
				end
			end
			
			ndir = this.lsf(l{:});
			ndir = reshape(ndir, Sign.count, Axis.count);
			ndir = -diff(ndir) ./ dl;			
			
			norm_n = norm(ndir);
			if norm_n < 10*eps
				ndir = [0 0 0];
			else
				ndir = ndir ./ norm_n;
			end
		end
		
% 		function [rvol, ndir] = smoothing_params4(this, box)
% 			chkarg(istypesizeof(box, 'real', [Axis.count, Sign.count]), ...
% 				'"box" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
% 
% 			l_probe = cell(1, Axis.count);	
% 			dl_sub = diff(box.') ./ this.n_subpxls;
% 			N_subpxls = this.n_subpxls^3;
% 			
% 			for w = Axis.elems
% 				bnd = box(w, :);  % interval in w-direction
% 
% % 				% Below, probing points are at the centers of subpixels.  One
% % 				% extra subpixel beyond each boundary of the box is considered
% % 				% to discard single-sided differences in gradient().
% % 				l_probe{w} = linspace(bnd(Sign.n) - dl_sub(w)/2, bnd(Sign.p) + dl_sub(w)/2, this.n_subpxls + 2);
% 
% 				l_probe{w} = linspace(bnd(Sign.n) + dl_sub(w)/2, bnd(Sign.p) - dl_sub(w)/2, this.n_subpxls);
% 			end
% 			
% 			[X, Y, Z] = ndgrid(l_probe{:});
% 			F = this.lsf(X, Y, Z);
% 			
% % 			Fin = F(2:end-1, 2:end-1, 2:end-1);  % values at interior points of box
% % 			is_contained = (Fin > 0);
% % 			at_interface = (Fin == 0);
% 
% 			is_contained = (F > 0);
% 			at_interface = (F == 0);
% 			rvol = sum(double(is_contained(:))) + 0.5 * sum(double(at_interface(:)));
% 			rvol = rvol / N_subpxls;
% 			
% % 			[Fy, Fx, Fz] = gradient(F, dl_sub(Axis.y), dl_sub(Axis.x), dl_sub(Axis.z));  % x and y are swapped due to MATLAB's definiton of gradient
% % 			Fx = Fx(2:end-1, 2:end-1, 2:end-1);
% % 			Fy = Fy(2:end-1, 2:end-1, 2:end-1);
% % 			Fz = Fz(2:end-1, 2:end-1, 2:end-1);
% % 			
% % 			gradF = [Fx(:), Fy(:), Fz(:)];
% % 			normF = sqrt(sum(gradF.^2, 2));
% % 			gradF = bsxfun(@rdivide, gradF, normF);
% % 			
% % 			ndir = sum(gradF) ./ N_subpxls;
% 
% 			dl = diff(box.');
% 			l = mat2cell(box.', Sign.count, ones(1, Axis.count));
% 			[X, Y, Z] = ndgrid(l{:});
% 			F = this.lsf(X, Y, Z);
% 			
% 			ndir = NaN(1, Axis.count);
% 			for w = Axis.elems
% 				dFw = -diff(F, 1, int(w));
% 				ndir(w) = mean(dFw(:)) / dl(w);
% 			end
% 			
% 			norm_n = norm(ndir);
% 
% 			% Each row of gradF has norm == 1.  If the norm of the average
% 			% gradient is much smaller than the norm of each gradient (== 1),
% 			% assume that the average gradient is zero; otherwise the error is
% 			% amplified by the division by a small number (= norm_n).
% 			if norm_n < 10*eps
% 				ndir = [0 0 0];
% 			else
% 				ndir = ndir ./ norm_n;
% 			end
% 		end
% 		
% 		function [rvol, ndir] = smoothing_params3(this, box)
% 			chkarg(istypesizeof(box, 'real', [Axis.count, Sign.count]), ...
% 				'"box" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
% 
% 			l_probe = cell(1, Axis.count);	
% 			dl_sub = diff(box.') ./ this.n_subpxls;
% 			N_subpxls = this.n_subpxls^3;
% 			
% 			for w = Axis.elems
% 				bnd = box(w, :);  % interval in w-direction
% 
% % 				% Below, probing points are at the centers of subpixels.  One
% % 				% extra subpixel beyond each boundary of the box is considered
% % 				% to discard single-sided differences in gradient().
% % 				l_probe{w} = linspace(bnd(Sign.n) - dl_sub(w)/2, bnd(Sign.p) + dl_sub(w)/2, this.n_subpxls + 2);
% 
% 				l_probe{w} = linspace(bnd(Sign.n) + dl_sub(w)/2, bnd(Sign.p) - dl_sub(w)/2, this.n_subpxls);
% 			end
% 			
% 			[X, Y, Z] = ndgrid(l_probe{:});
% 			F = this.lsf(X, Y, Z);
% 			
% % 			Fin = F(2:end-1, 2:end-1, 2:end-1);  % values at interior points of box
% % 			is_contained = (Fin > 0);
% % 			at_interface = (Fin == 0);
% 
% 			is_contained = (F > 0);
% 			at_interface = (F == 0);
% 			rvol = sum(double(is_contained(:))) + 0.5 * sum(double(at_interface(:)));
% 			rvol = rvol / N_subpxls;
% 			
% % 			[Fy, Fx, Fz] = gradient(F, dl_sub(Axis.y), dl_sub(Axis.x), dl_sub(Axis.z));  % x and y are swapped due to MATLAB's definiton of gradient
% % 			Fx = Fx(2:end-1, 2:end-1, 2:end-1);
% % 			Fy = Fy(2:end-1, 2:end-1, 2:end-1);
% % 			Fz = Fz(2:end-1, 2:end-1, 2:end-1);
% % 			
% % 			gradF = [Fx(:), Fy(:), Fz(:)];
% % 			normF = sqrt(sum(gradF.^2, 2));
% % 			gradF = bsxfun(@rdivide, gradF, normF);
% % 			
% % 			ndir = sum(gradF) ./ N_subpxls;
% 
% 			c = mean(box.');
% 			dl = diff(box.');
% 			l = mat2cell(repmat(c, Sign.count*Axis.count, 1), Sign.count*Axis.count, ones(1, Axis.count));
% 			for w = Axis.elems
% 				for s = Sign.elems
% 					l{w}(Sign.count * subsindex(w) + int(s)) = box(w, s);
% 				end
% 			end
% 			
% 			
% 			ndir = this.lsf(l{:});
% 			ndir = reshape(ndir, Sign.count, Axis.count);
% 			ndir = -diff(ndir) ./ dl;			
% 			
% 			norm_n = norm(ndir);
% 
% 			% Each row of gradF has norm == 1.  If the norm of the average
% 			% gradient is much smaller than the norm of each gradient (== 1),
% 			% assume that the average gradient is zero; otherwise the error is
% 			% amplified by the division by a small number (= norm_n).
% 			if norm_n < 10*eps
% 				ndir = [0 0 0];
% 			else
% 				ndir = ndir ./ norm_n;
% 			end
% 		end
% 		
% 		function [rvol, ndir] = smoothing_params1(this, box)
% 			chkarg(istypesizeof(box, 'real', [Axis.count, Sign.count]), ...
% 				'"box" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
% 
% 			l_probe = cell(1, Axis.count);	
% 			dl_sub = diff(box.') ./ this.n_subpxls;
% 			N_subpxls = this.n_subpxls^3;
% 			
% 			for w = Axis.elems
% 				bnd = box(w, :);  % interval in w-direction
% 
% 				% Below, probing points are at the centers of subpixels.  One
% 				% extra subpixel beyond each boundary of the box is considered
% 				% to discard single-sided differences in gradient().
% 				l_probe{w} = linspace(bnd(Sign.n) - dl_sub(w)/2, bnd(Sign.p) + dl_sub(w)/2, this.n_subpxls + 2);
% 			end
% 			
% 			[X, Y, Z] = ndgrid(l_probe{:});
% 			F = this.lsf(X, Y, Z);
% 			
% 			Fin = F(2:end-1, 2:end-1, 2:end-1);  % values at interior points of box
% 			is_contained = (Fin > 0);
% 			at_interface = (Fin == 0);
% 			rvol = sum(double(is_contained(:))) + 0.5 * sum(double(at_interface(:)));
% 			rvol = rvol / N_subpxls;
% 			
% 			[Fy, Fx, Fz] = gradient(F, dl_sub(Axis.y), dl_sub(Axis.x), dl_sub(Axis.z));  % x and y are swapped due to MATLAB's definiton of gradient
% 			Fx = Fx(2:end-1, 2:end-1, 2:end-1);
% 			Fy = Fy(2:end-1, 2:end-1, 2:end-1);
% 			Fz = Fz(2:end-1, 2:end-1, 2:end-1);
% 			
% 			gradF = [Fx(:), Fy(:), Fz(:)];
% 			normF = sqrt(sum(gradF.^2, 2));
% 			gradF = bsxfun(@rdivide, gradF, normF);
% 			
% 			ndir = -sum(gradF) ./ N_subpxls;  % (-) sign is for outward direction
% 
% 			norm_n = norm(ndir);
% 
% 			% Each row of gradF has norm == 1.  If the norm of the average
% 			% gradient is much smaller than the norm of each gradient (== 1),
% 			% assume that the average gradient is zero; otherwise the error is
% 			% amplified by the division by a small number (= norm_n).
% 			if norm_n < 10*eps
% 				ndir = [0 0 0];
% 			else
% 				ndir = ndir ./ norm_n;
% 			end
% 		end
% 		
% 		function [rvol, ndir] = smoothing_params2(this, box, center)
% 			chkarg(istypesizeof(box, 'real', [Axis.count, Sign.count]), ...
% 				'"box" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
% 			chkarg(istypesizeof(center, 'real', [1 Axis.count]), ...
% 				'"center" should be length-%d row vector with real elements.', Axis.count);
% 			chkarg(all(box(:,Sign.n) < center(:)) && all(box(:,Sign.p) > center(:)), ...
% 				'"center" should be included in "box".');
% 
% 			l_probe = cell(1, Axis.count);	
% 			dl_sub = diff(box.') ./ this.n_subpxls;
% 			N_subpxls = this.n_subpxls^3;
% 			
% 			for w = Axis.elems
% 				bnd = box(w, :);  % interval in w-direction
% 
% % 				% Below, probing points are at the centers of subpixels.  One
% % 				% extra subpixel beyond each boundary of the box is considered
% % 				% to discard single-sided differences in gradient().
% % 				l_probe{w} = linspace(bnd(Sign.n) - dl_sub(w)/2, bnd(Sign.p) + dl_sub(w)/2, this.n_subpxls + 2);
% 
% 				l_probe{w} = linspace(bnd(Sign.n) + dl_sub(w)/2, bnd(Sign.p) - dl_sub(w)/2, this.n_subpxls);
% 			end
% 			
% 			[X, Y, Z] = ndgrid(l_probe{:});
% 			F = this.lsf(X, Y, Z);
% 			
% % 			Fin = F(2:end-1, 2:end-1, 2:end-1);  % values at interior points of box
% % 			is_contained = (Fin > 0);
% % 			at_interface = (Fin == 0);
% 
% 			is_contained = (F > 0);
% 			at_interface = (F == 0);
% 			rvol = sum(double(is_contained(:))) + 0.5 * sum(double(at_interface(:)));
% 			rvol = rvol / N_subpxls;
% 			
% % 			[Fy, Fx, Fz] = gradient(F, dl_sub(Axis.y), dl_sub(Axis.x), dl_sub(Axis.z));  % x and y are swapped due to MATLAB's definiton of gradient
% % 			Fx = Fx(2:end-1, 2:end-1, 2:end-1);
% % 			Fy = Fy(2:end-1, 2:end-1, 2:end-1);
% % 			Fz = Fz(2:end-1, 2:end-1, 2:end-1);
% % 			
% % 			gradF = [Fx(:), Fy(:), Fz(:)];
% % 			normF = sqrt(sum(gradF.^2, 2));
% % 			gradF = bsxfun(@rdivide, gradF, normF);
% % 			
% % 			ndir = sum(gradF) ./ N_subpxls;
% 
% 			vc = [center - dl_sub/2; center + dl_sub/2];  % vicinity of center
% 			lc = mat2cell([center; center], 2, [1 1 1]);
% 			
% 			ndir = NaN(1, Axis.count);
% 			for w = Axis.elems
% 				loc = lc;
% 				loc{w} = vc(:, w);
% 				ndir(w) = -diff(this.lsf(loc{:})) ./ dl_sub(w);  % (-) sign is for outward direction
% 			end
% 			
% 			norm_n = norm(ndir);
% 
% 			% Each row of gradF has norm == 1.  If the norm of the average
% 			% gradient is much smaller than the norm of each gradient (== 1),
% 			% assume that the average gradient is zero; otherwise the error is
% 			% amplified by the division by a small number (= norm_n).
% 			if norm_n < 10*eps
% 				ndir = [0 0 0];
% 			else
% 				ndir = ndir ./ norm_n;
% 			end
% 		end		
	end
end
