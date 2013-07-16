% Sign is the enumeration class for negative and positive signs.
classdef Sign < Enumerated	
	enumeration
		n('negative')
		p('positive')
	end
	
	methods (Static)
		function elems = elems(ind)
			elems = [Sign.n, Sign.p];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(Sign.elems);
		end
	end
	
	methods
		function g = alter(this)
			% This function simply returned Sign.p if this == Sign.n and Sign.n
			% if this == Sign.p originally, but it is modified to handle
			% arguments given as arrays.
			chkarg(istypeof(this, 'Sign'), '"this" should have Sign as elements.');
			if isempty(this)  % alter([]) makes no sense, but alter(Sign.empty()) does
				g = [];
			else
				g = this.elems(GT.count + 1 - int(this));
			end
		end
	end
end

