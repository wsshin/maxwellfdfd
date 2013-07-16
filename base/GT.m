% GT is the enumeration class for the types of grids.
classdef GT < Enumerated
	enumeration
		prim('primary')
		dual('dual')
	end
	
	methods (Static)
		function elems = elems(ind)
			elems = [GT.prim, GT.dual];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(GT.elems);
		end
	end

	methods
		function g = alter(this)
			% This function simply returned GT.dual if this == GT.prim and
			% GT.prim if this == GT.dual originally, but it is modified to
			% handle arguments given as arrays.
			chkarg(istypeof(this, 'GT'), '"this" should have GT as elements.');
			if isempty(this)  % alter([]) makes no sense, but alter(GT.empty()) does
				g = [];
			else
				g = this.elems(GT.count + 1 - int(this));
			end
		end
	end
end

