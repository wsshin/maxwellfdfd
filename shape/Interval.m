classdef Interval
	% Interval is a 1D version of Shape.  It is a concrete class.
	
	properties (SetAccess = immutable)
		lprim  % lprim_array that is used to generate grid
		bound  % [min max];
		dl_max  % maximum dl
	end
	
	methods
		function this = Interval(lprim_array, dl_max)
			% lprim
			chkarg(istypesizeof(lprim_array, 'real', [1 0]) && ~isempty(lprim_array), ...
				'"lprim_array" should be nonempty row vector with real elements.');
			this.lprim = unique(lprim_array);  % sorted and duplicate elements are removed
			
			% bounds
			this.bound = [min(lprim_array), max(lprim_array)];
						
			% dl_max
			chkarg(istypesizeof(dl_max, 'real') && dl_max > 0, '"dl_max" should be postive and real.');
			this.dl_max = dl_max;
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
