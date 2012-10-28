% BC is the enumeration class for boundary conditions.
classdef BC < Enumerated
	enumeration
		e('PEC')  % tangential component of E-field = 0
		m('PMC')  % normal component of E-field = 0
		p('periodic')  % periodic boundary condition
	end

	methods (Static)
		function elems = elems(ind)
			elems = [BC.e, BC.m, BC.p];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(BC.elems);
		end
	end
end
