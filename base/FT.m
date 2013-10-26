% FT is the enumeration class for the typel of fields.
classdef FT < Enumerated
	enumeration
		e('E')
		h('H')
	end
	
	methods (Static)
		function elems = elems(ind)
			elems = [FT.e, FT.h];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(FT.elems);
		end
	end

	methods
		function f = alter(this)
			% This function simply returned FT.h if this == FT.e and FT.e if
			% this == FT.h originally, but it is modified to handle arguments
			% given as arrays.
			chkarg(istypeof(this, 'FT'), '"this" should have FT as elements.');
			if isempty(this)  % alter([]) makes no sense, but alter(FT.empty()) does
				f = [];
			else
				f = this.elems(FT.count + 1 - int(this));
			end
		end
	end
end

