classdef Material
	% Material is a class for a material.
	
	properties (SetAccess = immutable)
		name
		color
		eps
		mu
	end
	
	methods
		function this = Material(name, color, eps, mu)
			chkarg(istypesizeof(name, 'char', [1 0]), '"name" should be string.');
			this.name = name;
			
			chkarg(istypeof(color, 'char') ...
				|| (istypesizeof(color, 'real', [1 3]) && all(color <= 1) && all(color >= 0)), ...
				'"color" should be string or [r g b].');
			this.color = color;
			
			chkarg(istypesizeof(eps, 'complex'), '"eps" should be complex.');
			this.eps = eps;
			
			if nargin < 4  % no mu
				mu = 1.0;
			end
			chkarg(istypesizeof(mu, 'complex'), '"mu" should be complex.');
			this.mu = mu;
		end
	end
	
	methods (Static)
		function material = create(name, color, osc)
			chkarg(istypesizeof(name, 'char', [1 0]), '"name" should be string.');
			chkarg(istypesizeof(color, 'char', [1 0]) ...
				|| (istypesizeof(color, 'real', [1 3]) && all(color <= 1) && all(color >= 0)), ...
				'"color" should be string or [r g b].');
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			eV = osc.in_eV();
			
			maxwell_root = fileparts(fileparts(mfilename('fullpath')));
			param_dir = [maxwell_root, filesep, 'material', filesep];
			param_file = [name, '.mat'];
			param = load([param_dir, param_file]);  % eV, n, k are loaded
			chkarg(eV >= param.eV(1) && eV <= param.eV(end), '"eV" should be in the range described by %s', param_file);
			
			n = interp1(param.eV, param.n, eV);
			k = interp1(param.eV, param.k, eV);
			epsilon = (n - 1i * k)^2;
			material =  Material(name, color, epsilon);
		end
	end
end

