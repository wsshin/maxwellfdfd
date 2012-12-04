%% Material
% Class representing an electromagnetic material.

%%% Description
% |Material| represents an electromagnetic material with electric permittivity
% and magnetic permeability.  Users should provide the name and color of the
% material in the constructor

%%% Construction
%  mat = Material(name, color, eps)
%  mat = Material(name, color, eps, mu)
% 
% *Input Arguments*
%
% * |name|: name of the material.  Used to distinguish materials and label them
% in figures.
% * |color|: color of the material.  Used to visualize the material in figures.
% It can be MATLAB's reserved color code characters (e.g., |'r'| for red, |'k'|
% for black) or an RGB code |[r g b]| where |r|, |g|, |b| are real numbers
% between 0 and 1.  Use |color = 'none'| to hide the objects made of this
% material in figures.
% * |eps|: electric permittivity of the material.
% * |mu|: magnetic permeability of the material.  If not assigned, the default
% value |mu = 1| is used.

%%% Example
%   % Create an instance of Material.
%   vacuum =  Material('vacuum', 'none', 1.0);
%
%   % Use the constructed material in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', vacuum, Box([0 100; 0 50; 0 20]), {REMAINING ARGUMENTS});

%%% See Also
% <maxwell_run.html |maxwell_run|>

classdef Material
	
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
	
	methods
		function [sorted, ind] = sort(this, varargin)
			% varargin is the optional parameters (such as 'descend') of sort().
			n = length(this);
			names = cell(1, n);
			for i = 1:n
				names{i} = this(i).name;
			end
			[~, ind] = sort(names, varargin{:}); 
			sorted = this(ind);
		end
		
		function truth = ne(this, another)
			chkarg(all(size(this) == size(another)), '"this" and "another" should have same size.');
			dims = size(this);
			n = numel(this);
			this = this(:);
			another = another(:);
			truth = true(n,1);
			for i = 1:n
				truth(i)= ~isequal(this(i).name, another(i).name);
			end
			truth = reshape(truth, dims);
		end
	end
end

