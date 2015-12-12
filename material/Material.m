%% Material
% Class representing an electromagnetic material.

%%% Description
% |Material| represents an electromagnetic material with electric permittivity
% and magnetic permeability.  Users should provide the name and color of the
% material in the constructor

%%% Construction
%  mat = Material(name, color, eps, [islossless])
%  mat = Material(name, color, eps, mu, [islossless])
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
% * |mu|: magnetic permeability of the material.  If unassigned, the default
% value |mu = 1| is used.
% * |islossless|: optional argument.  |true| to ignore the loss or gain in
% material; |false| otherwise.  The default value is |false|.

%%% Example
%   % Create an instance of Material.
%   vacuum =  Material('vacuum', 'none', 1.0);
%
%   % Use the constructed material in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', vacuum, Box([0 100; 0 50; 0 20]), {REMAINING ARGUMENTS});

%%% See Also
% <maxwell_run.html |maxwell_run|>

classdef Material < handle
	
	properties (SetAccess = immutable)
		eps
		mu
	end

	properties (Dependent, SetAccess = immutable)
		hasisoeps
		hasisomu
		isiso
	end
	
	properties
		name
		color
	end

	methods
		function this = Material(name, color, eps, varargin)
			narginchk(3, 5);
			chkarg(istypesizeof(name, 'char', [1 0]), '"name" should be string.');
			chkarg(istypeof(color, 'char') ...
				|| (istypesizeof(color, 'real', [1 3]) && all(color <= 1) && all(color >= 0)), ...
				'"color" should be string or [r g b].');
			chkarg(istypeof(eps, 'complex') ...
				&& (isexpandable2row(eps, Axis.count) || isequal(size(eps), [Axis.count Axis.count]) ...
				|| isequal(size(eps), [Axis.count 2*Axis.count])), ...
				'"eps" should be complex scalar, length-%d row vector, %d-by-%d matrix, or %d-by-%d matrix.', ...
				Axis.count, Axis.count, Axis.count, Axis.count, 2*Axis.count);
			
			% Make eps a 3-by-3 tensor
			if isexpandable2row(eps, Axis.count)
				eps = diag(expand2row(eps, Axis.count));
			elseif isequal(size(eps), [Axis.count 2*Axis.count])
				S = eps(:, (Axis.count+1):2*Axis.count);  % columns of S are crystal axes
				eps = eps(:, 1:Axis.count);
				eps = S * eps / S;
			end
			
			mu_temp = 1.0;
			islossless = false;
			if length(varargin) == 2
				mu_temp = varargin{1};
				islossless = varargin{2};
			elseif length(varargin) == 1
				if istypesizeof(varargin{1}, 'logical')
					islossless = varargin{1};
				else
					mu_temp = varargin{1};
				end
			end
			
			chkarg(istypeof(mu_temp, 'complex') ...
				&& (isexpandable2row(mu_temp, Axis.count) || isequal(size(mu_temp), [Axis.count Axis.count]) ...
				|| isequal(size(mu_temp), [Axis.count 2*Axis.count])), ...
				'"mu" should be complex scalar, length-%d row vector, %d-by-%d matrix, or %d-by-%d matrix.', ...
				Axis.count, Axis.count, Axis.count, Axis.count, 2*Axis.count);

			% Make mu a 3-by-3 tensor
			if isexpandable2row(mu_temp, Axis.count)
				mu_temp = diag(expand2row(mu_temp, Axis.count));
			elseif isequal(size(mu_temp), [Axis.count 2*Axis.count])
				S = mu_temp(:, (Axis.count+1):2*Axis.count);  % columns of S are crystal axes
				mu_temp = mu_temp(:, 1:Axis.count);
				mu_temp = (S * mu_temp) / S;
			end
			
			chkarg(istypesizeof(islossless, 'logical'), '"islossless" should be logical.');
			
			if islossless
				name = [name, ' (lossless)'];
				eps = real(eps);
				mu_temp = real(mu_temp);
			end
			this.name = name;
			this.color = color;
			this.eps = eps;
			this.mu = mu_temp;
		end
		
		function truth = get.hasisoeps(this)
			truth = isdiag(this.eps) && all(diag(this.eps) == this.eps(1,1));
		end
		
		function truth = get.hasisomu(this)
			truth = isdiag(this.eps) && all(diag(this.mu) == this.mu(1,1));
		end
		
		function truth = get.isiso(this)
			truth = this.hasisoeps && this.hasisomu;
		end
	end
	
	methods (Static)
		function material = fromtable(osc, name, color, islossless)
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			chkarg(istypesizeof(name, 'char', [1 0]), '"name" should be string.');
			chkarg(istypesizeof(color, 'char', [1 0]) ...
				|| (istypesizeof(color, 'real', [1 3]) && all(color <= 1) && all(color >= 0)), ...
				'"color" should be string or [r g b].');
			
			if nargin < 4
				islossless = false;
			end
			chkarg(istypesizeof(islossless, 'logical'), '"islossless" should be logical.');
			
			eV = osc.in_eV();
			
			maxwell_root = fileparts(fileparts(mfilename('fullpath')));
			param_dir = [maxwell_root, filesep, 'dielconst', filesep];
			param_file = [strrep(name, '/', filesep), '.mat'];
			param = load([param_dir, param_file]);  % eV, n, k are loaded
			
			if eV < param.eV(1) && eV * (1+1e-6) >= param.eV(1)
				eV = param.eV(1);
			elseif eV > param.eV(end) && eV * (1-1e-6) <= param.eV(end)
				eV = param.eV(end);
			end
			chkarg(eV >= param.eV(1) && eV <= param.eV(end), '"eV" should be in the range described by %s', param_file);
			
			% Interpolate n and k rather than eps' and eps", because n and k are
			% measured data.
			n = interp1(param.eV, param.n, eV);
			k = interp1(param.eV, param.k, eV);
			epsilon = (n - 1i * k)^2;
			material =  Material(name, color, epsilon, islossless);
		end
	end

	% Need the following definitions to eliminate duplication in 'materials
	% used' information printed out in build_system().
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
				truth(i)= ~isequal(this(i).name, another(i).name) || ...
						~isequal(this(i).color, another(i).color) || ...
						any(this(i).eps(:) ~= another(i).eps(:)) || ...
						any(this(i).mu(:) ~= another(i).mu(:));
			end
			truth = reshape(truth, dims);
		end
	end
end

