% PML is the enumeration class for the kinds of PML
classdef PML < Enumerated	
	enumeration
		SC('stretched-coordinate')
		U('uniaxial')
	end

	methods (Static)
		function elems = elems(ind)
			elems = [PML.SC, PML.U];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(PML.elems);
		end
	end
end
