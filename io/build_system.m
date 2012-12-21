%% build_system
% Run MaxwellFDS.

%%% Syntax
%  [osc, grid3d, s_factor_cell, eps_face_cell, mu_edge_cell, J_cell] = build_system(OSC, DOM, OBJ, SRC, [progmark])
%  [..., obj_array, src_array] = build_system(OSC, DOM, OBJ, SRC, [pragmark])
%  [..., eps_node_array, mu_node_array] = build_system(OSC, DOM, OBJ, SRC, [pragmark])


%%% Description
% |build_system(OSC, DOM, OBJ, SRC, [progmark])| constructs a system from
% given objects and sources.  The constructed system is typically used inside
% <maxwell_run.html maxwell_run>.
%
% Each of |OSC|, |DOM|, |OBJ|, and |SRC| represents a group of parameters.
% Each group supports several flexible expressions.  For more details, see the
% relevant sections about the input parameter groups in <maxwell_run.html
% |maxwell_run|>.
%
% An additional input parameter |progmark| is an instance of <ProgMark.html
% ProgMark>, which outputs the progress of the system build procedure as the
% standard output.  If it is not given, then it is created internally.
%
% |[osc, grid3d, s_factor_cell, eps_face_cell, mu_edge_cell, J_cell] =
% build_system(...)| returns
%
% * |osc|, an instance of <Oscillation.html Oscillation>
% * |grid3d|, an instance of <Grid3d.html Grid3d>, 
% * |s_factor_cell|, a cell array of PML s-factors: |{sx_array, sy_array,
% sz_array}|
% * |eps_face_cell|, a cell array of electric permittivity evaluated at the
% centers of the finite-difference grid faces: |{eps_xx_array, eps_yy_array,
% eps_zz_array}|
% * |mu_edge_cell|,  a cell array of magnetic permeability evaluated at the
% centers of the finite-difference grid edges: |{mu_xx_array, mu_yy_array,
% mu_zz_array}|
% * |J_cell|, a cell array of electric current sources: |{Jx_array, Jy_array,
% Jz_array}|
% 
% |[..., obj_array, src_array] = build_system(...)| returns additionally arrays
% of instances of <Object.html |Object|> and <Source.html |Source|>.  The
% |Object| and |Source| elements represent the objects and sources placed in the
% simulation domain, so they can be used to visualize the simulation domain.
%
% |[..., eps_node, mu_node] = build_system(...)| returns additionally arrays of
% electric permittivity and magnetic permeability evaluated at the nodes of the
% finite-difference grid.


%%% Example
%   gray = [0.5 0.5 0.5];  % [r g b]
%   inspect_only = true;
%   [E, H, obj_array, err] = build_system(...
%       'OSC', 1e-9, 1550, ...
%       'DOM', {['Palik', filesep, 'SiO2'], 'none'}, [-700, 700; -600, 600; -200, 1700], 20, BC.p, 200, ...
%       'OBJ', ...
%           {['Palik', filesep, 'SiO2'], 'none'}, Box([-50, 50; -50, 50; -200, 1700], [2, 2, 20]), ...  % OBJ1
%           {['CRC', filesep, 'Ag'], gray}, [Box([-700, -25; -25, 25; -200, 1700], 20), Box([25, 700; -25, 25; -200, 1700], 20)], ...  % OBJ2
%       'SRC', PointSrc(Axis.x, [0, 0, 200]) ...
%       );

function [osc, grid3d, s_factor_cell, eps_face_cell, mu_edge_cell, J_cell, ...
	obj_array, src_array, eps_node_array, mu_node_array] = build_system(varargin)
	iarg = nargin; arg = varargin{iarg};
	if istypesizeof(arg, 'ProgMark')
		pm = arg;
		narglim = nargin - 1;
	else
	 	pm = ProgMark();
		narglim = nargin;
	end
		
	function material = create_material(varargin)
		if istypesizeof(varargin{end}, 'logical')
			islossless = varargin{end};
			varargin = varargin(1:end-1);
		else
			islossless = false;
		end
		
		if length(varargin) == 2
			material = Material.create(varargin{:}, osc, islossless);
		else
			material = Material(varargin{:});
		end
	end

	osc = Oscillation.empty();
	obj_dom = Object.empty();
	shape_array = Shape.empty();
	sshape_array = Shape.empty();
	mat_array = Material.empty();
	obj_array = Object.empty();
	sobj_array = Object.empty();
	src_array = [];
	isepsgiven = false;
	isTFSF = false;
	iarg = 0;
	while iarg < narglim
		iarg = iarg + 1; arg = varargin{iarg};
		
		if ischar(arg) && strcmpi(arg,'OSC')
			% Set up OSC.
			iarg = iarg + 1; arg = varargin{iarg};
			chkarg((istypesizeof(arg, 'real') && arg > 0) || istypesizeof(arg, 'Oscillation'), ...
				'"argument #%d should be either "L0" (positive) or "osc" (instance of Oscillation).', iarg);
			if istypesizeof(arg, 'real')
				L0 = arg;
				iarg = iarg + 1; wvlen = varargin{iarg};
				chkarg(istypesizeof(wvlen, 'complex'), 'argument #%d should be "wvlen" (complex).', iarg);
				unit = PhysUnit(L0);
				osc = Oscillation(wvlen, unit);
			else  % arg is instance of Oscillation
				osc = arg;
			end
		elseif ischar(arg) && strcmpi(arg,'DOM')
			% Set up DOM.
			iarg = iarg + 1; arg = varargin{iarg};
			chkarg(iscell(arg) || istypesizeof(arg, 'Material') || istypesizeof(arg, 'Object'), ...
				'argument #%d should be cell, instance of Material, or instance of Object.', iarg);
			if istypesizeof(arg, 'Object')
				obj_dom = arg;
				domain = obj_dom.shape;
				chkarg(istypesizeof(domain, 'Domain'), 'argument #%d should be instance of Object with Domain as its shape.', iarg);
			else
				if iscell(arg)
					mat_dom = create_material(arg{:});
				else
					assert(istypesizeof(arg, 'Material'));
					mat_dom = arg;
				end
		
				iarg = iarg + 1; arg = varargin{iarg};
				chkarg(istypesizeof(arg, 'real', [Axis.count, Sign.count]) || istypesizeof(arg, 'Domain'), ...
					'argument #%d should be either "box_dom" ([xmin xmax; ymin ymax; zmin zmax]) or "domain" (instance of Domain).', iarg);
				if istypesizeof(arg, 'real', [Axis.count, Sign.count])
					box_domain = arg;
					iarg = iarg + 1; dl_domain = varargin{iarg};
					chkarg(istypeof(dl_domain, 'real') && isexpandable2row(dl_domain, Axis.count), ...
						'"argument #%d should be dl_domain (positive number or length-%d row vector of positive numbers).', iarg, Axis.count);
					domain = Domain(box_domain, expand2row(dl_domain, Axis.count));
				else  % arg is instance of Domain
					domain = arg;
				end
				obj_dom = Object(domain, mat_dom);
			end
			mat_array = [mat_array(1:end), obj_dom.material];
	
			% Set up boundary conditions and PML thicknesses.
			iarg = iarg + 1; bc = varargin{iarg};
			chkarg(istypeof(bc, 'BC') && isexpandable2mat(bc, Axis.count, Sign.count), ...
				'argument #%d should be "bc" (scalar, length-%d row vector, or %d-by-%d matrix with BC as elements).', iarg, Axis.count, Axis.count, Sign.count);
			iarg = iarg + 1; Lpml = varargin{iarg};
			chkarg(istypeof(Lpml, 'real') && isexpandable2mat(Lpml, Axis.count, Sign.count) && all(all(Lpml>=0)), ...
				'argument #%d should be "Lpml" (scalar, length-%d row vector, or %d-by-%d matrix with nonnegative numbers as elements).', iarg, Axis.count, Axis.count, Sign.count);
	
			% Set up a flag to generate a grid dynamically.
			iarg = iarg + 1; arg = varargin{iarg};
			withuniformgrid = false;  % generate a grid dynamically by default.
			if istypesizeof(arg, 'logical')
				withuniformgrid = arg;
			else
				iarg = iarg - 1; % because withuniformgrid is optional argument
			end
		elseif ischar(arg) && (strcmpi(arg,'OBJ') || strcmpi(arg,'SOBJ'))
			% Set up OBJ.
			is_scatterer = strcmpi(arg,'SOBJ');
			iarg = iarg + 1; arg = varargin{iarg};
			if istypesizeof(arg, 'complex', [0 0 0])
				isepsgiven = true;
				eps_node_array = arg;
				mu_node_array = ones(size(eps_node_array));
			else
				% Set up objects.
				obj_array_temp = Object.empty();
				shape_array_temp = Shape.empty();
				while iscell(arg) || istypesizeof(arg, 'Material') || istypesizeof(arg, 'Object', [1 0])
					if istypesizeof(arg, 'Object', [1 0])
						objs = arg;
						obj_array_temp = [obj_array_temp(1:end), objs];
						for obj = objs
							shape_array_temp = [shape_array_temp(1:end), obj.shape];
							mat_array = [mat_array(1:end), obj.material];
						end
						iarg = iarg + 1; arg = varargin{iarg};
					else
						if iscell(arg)
							mat = create_material(arg{:});
						else
							assert(istypesizeof(arg, 'Material'));
							mat = arg;
						end
						mat_array = [mat_array(1:end), mat];

						iarg = iarg + 1; arg = varargin{iarg};
						while istypesizeof(arg, 'Shape', [1 0])
							shapes = arg;
							shape_array_temp = [shape_array_temp(1:end), shapes];
							objs = Object(shapes, mat);
							obj_array_temp = [obj_array_temp(1:end), objs];
							iarg = iarg + 1; arg = varargin{iarg};
						end
					end
				end
				iarg = iarg - 1;
				
				if is_scatterer
					sshape_array = [sshape_array(1:end), shape_array_temp];
					sobj_array = [sobj_array(1:end), obj_array_temp];
				else
					shape_array = [shape_array(1:end), shape_array_temp];
					obj_array = [obj_array(1:end), obj_array_temp];
				end
			end
		elseif ischar(arg) && strcmpi(arg,'SRC')
			% Set up sources.
			iarg = iarg + 1; arg = varargin{iarg};
			if ~istypesizeof(arg, 'Source', [1 0])
				warning('FDS:buildSys', 'no source is given.');
			end

			while istypesizeof(arg, 'Source', [1 0])
				if istypesizeof(arg, 'TFSFPlaneSrc')
					isTFSF = true;
				end
				src_array = [src_array(1:end), arg];
				iarg = iarg + 1; arg = varargin{iarg};
			end
			iarg = iarg - 1;
		elseif iarg == narglim
			chkarg(false, ['some arguments are not used.\n', ...
				'Suggestion: check if each parameter group is specified with beginning specifier.']);
		end
	end
	chkarg(~isempty(osc), 'OSC parameter groups should be set.');
	chkarg(~isempty(obj_dom), 'DOM parameter groups should be set.');
	obj_array = [obj_dom, obj_array];
	
	chkarg(iarg <= narglim, 'more arguments than expected.');	
	pm.mark('initial setup');
	
	mat_array = unique(mat_array);
	fprintf('materials used:\n');
	for mat = mat_array
		fprintf('\t%s: eps = %s, mu = %s\n', mat.name, num2str(mat.eps), num2str(mat.mu));
	end

	% Generate a grid.
	[lprim, Npml] = generate_lprim3d(domain, Lpml, [shape_array, sshape_array], src_array, withuniformgrid);
	grid3d = Grid3d(osc.unit, lprim, Npml, bc);
	if withuniformgrid
		pm.mark('uniform grid generation');
	else
		pm.mark('dynamic grid generation');
	end
	fprintf('\t[Nx Ny Nz] = %s\n', mat2str(grid3d.N));


	% Construct material parameters.
	if ~isepsgiven
		[eps_node_array, mu_node_array] = assign_material_node(grid3d, obj_array);  % (Nx+1) x (Ny+1) x (Nz+1)
	end
	eps_face_cell = harmonic_mean_eps_node(eps_node_array);
	mu_edge_cell = arithmetic_mean_mu_node(mu_node_array);

	% Construct PML s-factors.
	s_factor_cell = generate_s_factor(osc.in_omega0(), grid3d);
	pm.mark('eps and mu assignment');

	if ~isTFSF
		% Solve for modes.
		for src = src_array
			if istypesizeof(src, 'ModalSrc')
				modalsrc = src;
				if ~modalsrc.ispreped
					prep_modalsrc(osc, grid3d, eps_face_cell, mu_edge_cell, s_factor_cell, modalsrc);
				end

				neff = modalsrc.neff;
				beta = 2*pi*neff / osc.in_L0();
				pm.mark('mode calculation');
				fprintf('\tbeta = %s, n_eff = %s\n', num2str(beta), num2str(neff));
			end
		end
	else  % isTFSF == true
		% Set up J for TF/SF.
		for src = src_array
			if istypesizeof(src, 'TFSFPlaneSrc')
				tfsfsrc = src;
				tfsfsrc.set_bg_material(obj_dom.material);
				E0 = tfsfsrc.create_incidentE(osc, grid3d);
				J = cell(1, Axis.count);
				for w = Axis.elems
					J{w} = zeros(grid3d.N);
				end
				
				d_prim = flip_vec(grid3d.dl(:, GK.dual));  % GK.dual, not GK.prim
				d_dual = flip_vec(grid3d.dl(:, GK.prim));  % GK.prim, not GK.dual
				s_prim = flip_vec(s_factor_cell(:, GK.dual));  % GK.dual, not GK.prim
				s_dual = flip_vec(s_factor_cell(:, GK.prim));  % GK.prim, not GK.dual
				mu_edge_temp = flip_vec(mu_edge_cell);
				eps_face_temp = flip_vec(eps_face_cell);
				J = neg_vec(flip_vec(J));  % pseudovector
				E0 = neg_vec(flip_vec(E0));  % pseudovector
	
				nosolve = true;
				[~, ~, A] = solve_eq_direct(osc.in_omega0(), ...
								d_prim, d_dual, ...
								s_prim, s_dual, ...
								mu_edge_temp, eps_face_temp, ...
								J, nosolve);
							
				x0 = [E0{Axis.x}(:); E0{Axis.y}(:); E0{Axis.z}(:)];
				r = reordering_indices(Axis.count, grid3d.N);
				x0 = x0(r);
							
				J = (A*x0) ./ (-1i*osc.in_omega0());
				J = reshape(J, [Axis.count grid3d.N]);
				J = permute(J, [int([Axis.x Axis.y Axis.z])+1, 1]);
				J = {J(:,:,:,Axis.x), J(:,:,:,Axis.y), J(:,:,:,Axis.z)};
				J = neg_vec(flip_vec(J));  % pseudovector
				
				tfsfsrc.setJ(J, grid3d);
				pm.mark('TF/SF source assignment');
			end
		end
		
		[eps_node_array, mu_node_array] = assign_material_node(grid3d, sobj_array, ...
			eps_node_array(1:end-1, 1:end-1, 1:end-1), mu_node_array(1:end-1, 1:end-1, 1:end-1));  % (Nx+1) x (Ny+1) x (Nz+1)
		eps_face_cell = harmonic_mean_eps_node(eps_node_array);
		mu_edge_cell = arithmetic_mean_mu_node(mu_node_array);
	end
	obj_array = [obj_array, sobj_array];

	% Construct sources.
	J_cell = assign_source(grid3d, src_array);
	pm.mark('J assignment');
end
