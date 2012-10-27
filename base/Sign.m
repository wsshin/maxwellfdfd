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
end

