% EquationType defines the type of the equation.
classdef EquationType
	properties (SetAccess = immutable)
		f  % field which this equation solves for; instance of FT
		ge  % grid type of the E-field grid; instance of GT
	end
	
	methods
        function this = EquationType(f, ge)
			chkarg(istypesizeof(f, 'FT'), '"f" should be instance of FT.');
			chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');
			this.f = f;
			this.ge = ge;
		end
	end
end
