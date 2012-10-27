% GK is the enumeration class for the kinds of grids.
classdef GK < Enumerated
	enumeration
		prim('primary')
		dual('dual')
	end
	
	methods (Static)
		function elems = elems(ind)
			elems = [GK.prim, GK.dual];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(GK.elems);
		end
	end

	methods
		function g = alter(this)
			% This function simply returned GK.dual if this == GK.prim and
			% GK.prim if this == GK.dual originally, but it is modified to
			% handle arguments given as arrays.
			chkarg(istypeof(this, 'GK'), '"this" should have GK as elements.');
			if isempty(this)  % alter([]) makes no sense, but alter(GK.empty()) does
				g = [];
			else
				g = this.elems(GK.count + 1 - int(this));
			end
		end
	end
end

