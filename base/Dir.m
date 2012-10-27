% Dir is the enumeration class for the horizontal and vertical directions.
classdef Dir < Enumerated
	enumeration
		h('horizontal')
		v('vertical')
	end

	methods (Static)
		function elems = elems(ind)
			elems = [Dir.h, Dir.v];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(Dir.elems);
		end
	end
	
	methods
		function ortho = alter(this)
			% This function simply returned Dir.v if this == Dir.h and Dir.h if
			% this == Dir.v originally, but it is modified to handle arguments
			% given as arrays.
			chkarg(istypeof(this, 'Dir'), '"this" should have Dir as elements.');
			if isempty(this)  % alter([]) makes no sense, but alter(Dir.empty()) does
				ortho = [];
			else
				ortho = this.elems(Dir.count + 1 - int(this));
			end
		end		
	end
end

