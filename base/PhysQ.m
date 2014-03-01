% PhysQ is the enumeration class representing physical quantities.
classdef PhysQ < Enumerated
	enumeration
		arbitrary('arbitrary', '', 'AU')  % AU: arbitrary unit
		L('length', 'L', 'm')
		omega('frequency', '\omega', 'rad/s')
		eps('permittivity', '\epsilon', 'F/m')
		mu('permeability', '\mu', 'H/m')
		E('E-field', 'E', 'V/m')
		H('H-field', 'H', 'A/m')
		J('electric current density', 'J', 'A/m^2')
		M('magnetic current density', 'M', 'V/m^2')
		I('electric current', 'I', 'A')
		V('electric voltage', 'V', 'V')
		S('Poynting vector', 'S', 'W/m^2')
		P('power', 'P', 'W')
		u('power desnsity', 'u', 'W/m^3')
		% Whenever new quantity is added, add it to elems below.
	end
	
	properties (SetAccess = immutable)
		symbol
		SIunit
	end
	
	methods
		function this = PhysQ(name, symbol, SIunit)
			chkarg(ischar(SIunit), '"SIunit" should be string.');
			chkarg(ischar(symbol), '"symbol" should be string');
			this = this@Enumerated(name);
			this.symbol = symbol;
			this.SIunit = SIunit;
		end
	end
	
	methods (Static)
		function elems = elems(ind)
			elems = [PhysQ.arbitrary, PhysQ.L, PhysQ.omega, PhysQ.eps, PhysQ.mu, PhysQ.E, PhysQ.H, PhysQ.J, PhysQ.M, PhysQ.I, PhysQ.V, PhysQ.S, PhysQ.P PhysQ.u];
			if nargin > 0  % ind
				elems = elems(ind);
			end
		end
		
		function count = count()
			count = length(PhysQ.elems);
		end
	end
end

