% Oscillation is a class containing the information on the operation frequency.
classdef Oscillation < handle
    properties (SetAccess = immutable, GetAccess = private)
        wvlen  % vacuum wavelength
	end

    properties (SetAccess = immutable)
		unit  % physical unit
    end
	
	methods
        function this = Oscillation(wvlen, unit)
			chkarg(istypesizeof(wvlen, 'complex'), '"wvlen" should be complex.');
            this.wvlen = wvlen;
			chkarg(istypesizeof(unit, 'PhysUnit'), '"unit" should be instance of PhysUnit.');
			this.unit = unit;
		end
		
		function wvlen = in_L0(this)
			wvlen = this.wvlen;
		end
		        
        function omega = in_omega0(this)
            omega = 2*pi / this.wvlen;  % omega = angular_freq / unit.omega0
        end
        
        function energy = in_eV(this)
            energy = PhysC.h * PhysC.c0 / (this.wvlen * this.unit.value(PhysQ.L));  %  E (in eV) = h c / lambda
        end
	end
	
end

