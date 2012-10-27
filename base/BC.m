% BC is the enumeration class for boundary conditions.
classdef BC < Enumerated
	enumeration
		Et0('Et = 0')  % tangential component of E-field = 0
		En0('En = 0')  % normal component of E-field = 0
		p('periodic')  % periodic boundary condition
	end

	methods (Static)
		function elems = elems(ind)
			elems = [BC.Et0, BC.En0, BC.p];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(BC.elems);
		end
	end
end
