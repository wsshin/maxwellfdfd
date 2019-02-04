classdef PhysUnit < handle
	% PhysUnit defines the units of the physical quantities used in MaxwellFDFD.
	
	properties (SetAccess = immutable, GetAccess = public)
		va = NaN(1,PhysQ.count);  % array of values of units
	end
			
	methods
		function this = PhysUnit(L0)
			this.va(PhysQ.arbitrary) = 1;
			this.va(PhysQ.L) = L0;  % length in m (L0-dependent)
			this.va(PhysQ.omega) = PhysC.c0 / this.va(PhysQ.L); % frequency in rad/s (L0-dependent)
			this.va(PhysQ.eps) = PhysC.eps0;  % permittivity in eps0
			this.va(PhysQ.mu) = PhysC.mu0;  % permeability in mu0
			this.va(PhysQ.E) = 1;  % E-field in V/m
			this.va(PhysQ.H) = this.va(PhysQ.E) / PhysC.eta0;  % H-field in A/m
			this.va(PhysQ.J) = this.va(	PhysQ.H) / this.va(PhysQ.L);  % electric current density in A/m^2 (L0-dependent)
			this.va(PhysQ.M) = this.va(PhysQ.E) / this.va(PhysQ.L);  % magnetic current density in A/m^2 (L0-dependent)
			this.va(PhysQ.I) = this.va(PhysQ.J) * this.va(PhysQ.L)^2;  % electric current in Amperes (L0-dependent)
			this.va(PhysQ.V) = this.va(PhysQ.E) * this.va(PhysQ.L);  % voltage in Volts (L0-dependent)
			this.va(PhysQ.S) = this.va(PhysQ.E) * this.va(PhysQ.H);  % Poynting vector in Watt/m^2
			this.va(PhysQ.P) = this.va(PhysQ.S) * this.va(PhysQ.L)^2;  % power in Watt (L0-dependent)
			this.va(PhysQ.u) = this.va(PhysQ.S) / this.va(PhysQ.L);  % power density in Watt/m^3 (L0-dependent)
		end

		function v0 = value(this, physQcell)
			chkarg(istypesizeof(physQcell, 'PhysQ') || isphysQcell(physQcell), ...
				'"physQcell" should be instance of PhysQ, or cell array {PhysQ, int; PhysQ, int; ...}.');
			if istypesizeof(physQcell, 'PhysQ')
				physQcell = {physQcell, 1};
			end
			
			n = size(physQcell, 1);
			v0 = 1;
			for i = 1:n
				physQ = physQcell{i,1};
				physDim = physQcell{i,2};
				v0 = v0 * this.va(physQ)^physDim;
			end
		end
		
% 		function v_normal = SI2normal(this, v_SI, physQ)
% 			chkarg(istypesizeof(physQ, 'PhysQ'), '"physQ" should be instance of PhysQ.');
% 			v_normal = v_SI ./ this.va(physQ);
% 		end
% 		
% 		function v_SI = normal2SI(this, v_normal, physQ)
% 			chkarg(istypesizeof(physQ, 'PhysQ'), '"physQ" should be instance of PhysQ.');
% 			v_SI = v_normal .* this.va(physQ);
% 		end
	end
end

