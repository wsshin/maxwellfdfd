classdef Interval
	% Interval is a 1D version of Shape.  It is a concrete class.
	
	properties (SetAccess = immutable)
		bound  % [min max];
		dl_max  % maximum dl
		dl_boundary  % dl at boundaries; [dl_n dl_p];
	end
	
	methods
		function this = Interval(bound, dl_max, dl_boundary)
			% bounds
			chkarg(istypesizeof(bound, 'real', [1, Sign.count]), '"bound" should be [min max].');
			chkarg(bound(Sign.n) <= bound(Sign.p), ...
				'lower bound should be smaller than upper bound of "bound".');
			this.bound = bound;
						
			% dl_max
			if nargin < 2  % no dl_max
				this.dl_max = NaN;
				this.dl_boundary = NaN(1, Sign.count);
			else  % dl_max
				chkarg(istypesizeof(dl_max, 'real') && dl_max > 0, '"dl_max" should be postive and real.');
				this.dl_max = dl_max;

				% dl_boundary
				if nargin < 3  % no dl_boundary
					dl_boundary = this.dl_max;
				end
				chkarg(istypesizeof(dl_boundary, 'real', [1 0]) && all(dl_boundary > 0), 'elements of "dl_boundary" should be positive.');
				chkarg(isexpandable2row(dl_boundary, Sign.count), ...
					'"dl_boundary" should be scalar or length-%d vector.', Sign.count);
				dl_boundary = expand2row(dl_boundary, Sign.count);
				for s = Sign.elems
					chkarg(dl_boundary(s) <= this.dl_max, ...
						'elements of "dl_boundary" should be smaller than "dl_max".');
				end
				this.dl_boundary = dl_boundary;
			end
		end
						
		function L = L(this)
			L = this.bound(Sign.p) - this.bound(Sign.n);
		end
						
		function [truth, distance] = contains(this, val)
			chkarg(istypesizeof(val, 'real', [0 1]), '"val" should be column vector with real elements.');
			bn = this.bound(Sign.n);
			bp = this.bound(Sign.p);
			
			truth = (val >= bn) & (val <= bp);
			
			if nargout >= 2  % distance
				distance = min(abs([val-bn, val-bp]), [], 2);
			end
		end
	end
end
