% BC is the enumeration class for boundary conditions.
classdef BC < Enumerated
	enumeration
		Ht0('Ht = 0')  % tangential component of H-field = 0
		Hn0('Hn = 0')  % normal component of H-field = 0
		p('periodic')  % periodic boundary condition
	end

	methods (Static)
		function elems = elems(ind)
			elems = [BC.Ht0, BC.Hn0, BC.p];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(BC.elems);
		end
	end
end
