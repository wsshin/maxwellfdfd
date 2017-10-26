%% build_system
% Build a system to solve by MaxwellFDFD.

%%% Syntax
%  [osc, grid3d, s_factor_cell, eps_cell, mu_cell, J_cell] = build_system(ge, OSC, DOM, OBJ, SRC, [progmark])
%  [..., obj_array, src_array, mat_array] = build_system(ge, OSC, DOM, OBJ, SRC, [pragmark])
%  [..., eps_node_array, mu_node_array] = build_system(ge, OSC, DOM, OBJ, SRC, [pragmark])


%%% Description
% |build_system(ge, OSC, DOM, OBJ, SRC, [progmark])| constructs a system from
% given objects and sources.  The constructed system is typically used inside
% <maxwell_run.html maxwell_run>.
%
% |ge| is an instance of |GT| and indicates the grid type of the _E_-field.
% Each following argument, |OSC|, |DOM|, |OBJ|, and |SRC|, represents a
% group of parameters. Each group supports several flexible expressions.
% For more details, see the relevant sections about the input parameter
% groups in <maxwell_run.html |maxwell_run|>.
%
% An additional input parameter |progmark| is an instance of <ProgMark.html
% ProgMark>, which outputs the progress of the system build procedure as the
% standard output.  If it is not given, then it is created internally.
%
% |[osc, grid3d, s_factor_cell, eps_cell, mu_cell, J_cell, M_cell] = build_system(...)|
% returns
%
% * |osc|, an instance of <Oscillation.html Oscillation>
% * |grid3d|, an instance of <Grid3d.html Grid3d>, 
% * |s_factor_cell|, a cell array of PML s-factors: |{sx_array, sy_array,
% sz_array}|
% * |eps_cell|, a cell array of electric permittivity evaluated at the E-field
% positions: |{eps_xx_array, eps_yy_array, eps_zz_array}|
% * |mu_cell|,  a cell array of magnetic permeability evaluated at the H-field
% positions: |{mu_xx_array, mu_yy_array, mu_zz_array}|
% * |J_cell|, a cell array of electric current sources: |{Jx_array, Jy_array,
% Jz_array}|
% * |M_cell|, a cell array of electric current sources: |{Mx_array, My_array,
% Mz_array}|
% 
% |[..., obj_array, src_array, mat_array] = build_system(...)| returns
% additionally arrays of instances of <EMObject.html |EMObject|>, <Source.html
% |Source|>, and <Material.html |Material|>.  The |EMObject| and |Source|
% elements represent the objects and sources placed in the simulation
% domain, so they can be used to visualize the simulation domain.
%
% |[..., eps_node_array, mu_node_array] = build_system(...)| returns
% additionally arrays of electric permittivity and magnetic permeability
% evaluated at the nodes of the finite-difference grid.


%%% Example
%   gray = [0.5 0.5 0.5];  % [r g b]
%   [osc, grid3d, s_factor_cell, eps_cell, mu_cell, J_cell, M_cell,	...
%       obj_array, src_array, mat_array, eps_node_array, mu_node_array] = build_system(...
%       'OSC', 1e-9, 1550, ...
%       'DOM', {['Palik/SiO2'], 'none'}, [-700, 700; -600, 600; -200, 1700], 20, BC.p, 200, ...
%       'OBJ', ...
%           {['Palik/SiO2'], 'none'}, Box([-50, 50; -50, 50; -200, 1700], [2, 2, 20]), ...  % OBJ1
%           {['CRC/Ag'], gray}, [Box([-700, -25; -25, 25; -200, 1700], 20), Box([25, 700; -25, 25; -200, 1700], 20)], ...  % OBJ2
%       'SRC', PointSrc(Axis.x, [0, 0, 200]) ...
%       );

function [osc, grid3d, s_factor_cell, eps_cell, mu_cell, J_cell, M_cell, ...
	obj_array, src_array, mat_array, eps_node, mu_node, isiso] = build_system(varargin)

	iarg = nargin; arg = varargin{iarg};
	if istypesizeof(arg, 'ProgMark')
		pm = arg;
		narglim = nargin - 1;
	else
	 	pm = ProgMark();
		varargin = [varargin, {pm}];
		narglim = nargin;
	end
	
	iarg = 1; arg = varargin{iarg};
	chkarg(istypesizeof(arg, 'GT'), 'argument #%d should be "ge" (GT).', iarg);
	ge = arg;

	iarg = iarg + 1; arg = varargin{iarg};
	chkarg(istypesizeof(arg, 'PML'), 'argument #%d should be "pml" (PML).', iarg);
	pml = arg;
	
	function material = create_material(varargin)
		narg = nargin;
		if istypesizeof(varargin{end}, 'logical')
			narg = narg - 1;
		end
		
		chkarg(narg >= 2, '# of arguments should be at least 2.')
		matname = varargin{1};
		if isempty(strfind(matname, '/'))  % data table is not specified
			material = Material(varargin{:});
		else  % data table is specified
			material = Material.fromtable(osc, varargin{:});
		end
	end

	osc = Oscillation.empty();
	obj_dom = EMObject.empty();
	shape_array = Shape.empty();
	sshape_array = Shape.empty();
	mat_array = Material.empty();
	obj_array = EMObject.empty();
	sobj_array = EMObject.empty();
	srcj_array = [];
	srcm_array = [];
	isepsgiven = false;
	isTFSF = false;
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
			chkarg(iscell(arg) || istypesizeof(arg, 'Material') || istypesizeof(arg, 'EMObject'), ...
				'argument #%d should be cell, instance of Material, or instance of EMObject.', iarg);
			if istypesizeof(arg, 'EMObject')
				obj_dom = arg;
				domain = obj_dom.shape;
				chkarg(istypesizeof(domain, 'Domain'), 'argument #%d should be instance of EMObject with Domain as its shape.', iarg);
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
				obj_dom = EMObject(domain, mat_dom);
			end
			mat_array = [mat_array(1:end), obj_dom.material];
	
			% Set up boundary conditions and PML thicknesses.
			iarg = iarg + 1; bc = varargin{iarg};
			chkarg(istypeof(bc, 'BC') && isexpandable2mat(bc, Axis.count, Sign.count), ...
				'argument #%d should be "bc" (scalar, length-%d row vector, or %d-by-%d matrix with BC as elements).', iarg, Axis.count, Axis.count, Sign.count);
			iarg = iarg + 1; Lpml = varargin{iarg};
			chkarg(istypeof(Lpml, 'real') && isexpandable2mat(Lpml, Axis.count, Sign.count) && all(all(Lpml>=0)), ...
				'argument #%d should be "Lpml" (scalar, length-%d row vector, or %d-by-%d matrix with nonnegative numbers as elements).', iarg, Axis.count, Axis.count, Sign.count);

			% Set up the degree of the polynomial grading of the PML scale
			% factors.
			iarg = iarg + 1; arg = varargin{iarg};
			deg_pml = 4;  % polynomial degree
			if istypeof(arg, 'real')
				deg_pml = arg;
			else
				iarg = iarg - 1; % because deg_pml is optional argument
			end
	
			% Set up the target reflection coefficient of the PML.
			iarg = iarg + 1; arg = varargin{iarg};
			R_pml = exp(-16);  % target reflection coefficient
			if istypeof(arg, 'real')
				R_pml = arg;
			else
				iarg = iarg - 1; % because R_pml is optional argument
			end
			
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
			if istypesizeof(arg, 'complex', [0 0 0])  % 3D complex array with arbitrary size
				isepsgiven = true;
				eps_node_cell = {arg, arg, arg};
				mu_node_temp = ones(size(arg));
				mu_node_cell = {mu_node_temp, mu_node_temp, mu_node_temp};
			else
				% Set up objects.
				obj_array_temp = EMObject.empty();
				shape_array_temp = Shape.empty();
				while iscell(arg) || istypesizeof(arg, 'Material') || istypesizeof(arg, 'EMObject', [1 0])
					if istypesizeof(arg, 'EMObject', [1 0])
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
							objs = EMObject(shapes, mat);
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
		elseif ischar(arg) && strcmpi(arg,'SRCJ')
			% Set up sources.
			iarg = iarg + 1; arg = varargin{iarg};
			if ~istypesizeof(arg, 'Source', [1 0])
				warning('Maxwell:buildSys', 'no source is given.');
			end

			while istypesizeof(arg, 'Source', [1 0])
				if istypesizeof(arg, 'TFSFPlaneSrc')
					isTFSF = true;
				end
				srcj_array_curr = arg;
				for src = srcj_array_curr
					src.set_gridtype(ge);
				end
				srcj_array = [srcj_array(1:end), srcj_array_curr];
				iarg = iarg + 1; arg = varargin{iarg};
			end
			iarg = iarg - 1;
		elseif ischar(arg) && strcmpi(arg,'SRCM')
			% Set up sources.
			iarg = iarg + 1; arg = varargin{iarg};
			if ~istypesizeof(arg, 'Source', [1 0])
				warning('Maxwell:buildSys', 'no source is given.');
			end

			while istypesizeof(arg, 'Source', [1 0])
				if istypesizeof(arg, 'TFSFPlaneSrc')
					isTFSF = true;
				end
				srcm_array_curr = arg;
				for src = srcm_array_curr
					src.set_gridtype(alter(ge));
				end
				srcm_array = [srcm_array(1:end), srcm_array_curr];
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
	src_array = [srcj_array, srcm_array];
	if isTFSF && isempty(sobj_array)
		warning('Maxwell:objAssign', 'TF/SF source is used, but scatteres are not defined in SOBJ group.');
	end
	
	chkarg(iarg <= narglim, 'more arguments than expected.');	
	pm.mark('initial setup');
	
	fprintf('\tLength Unit: %s m\n', osc.unit.value(PhysQ.L));
	fprintf('\twvlen = %s, freq = %s eV\n', num2str(osc.in_L0()), num2str(osc.in_eV()));
	
	mat_array = unique(mat_array);
	isiso = true;
	fprintf('materials used:\n');
	for mat = mat_array
		epstext = mat.eps;
		if length(unique(epstext)) == 1
			epstext = epstext(Axis.x);
		end
		
		mutext = mat.mu;
		if length(unique(mutext)) == 1
			mutext = mutext(Axis.x);
		end
		
		fprintf('\t%s: eps = %s, mu = %s\n', mat.name, mat2str(epstext), mat2str(mutext));
		
		isiso = isiso && mat.isiso;
	end

	% Generate a grid.
	[lprim, Npml] = generate_lprim3d(domain, Lpml, [shape_array, sshape_array], src_array, withuniformgrid);
	grid3d = Grid3d(osc.unit, lprim, Npml, bc);
	if withuniformgrid
		pm.mark('uniform grid generation');
	else
		pm.mark('nonuniform grid generation');
	end
	fprintf('\t[Nx Ny Nz] = %s\n', mat2str(grid3d.N));
	
	% Generate a warning when a seemingly 2D simulation is defined on a 3D grid.
	[like2d, normal_axis] = is2dlike(grid3d.N);
	if like2d && grid3d.N(normal_axis) >= 2  % possible user mistake
		warning(['If this a 2D structure, N%s should be 1; ', ...
			'check d%s''s of objects and locations of sources'], char(normal_axis), char(normal_axis));
	end

	% Construct material parameters.
	if ~isepsgiven
		[eps_node_cell, mu_node_cell] = assign_material_node(grid3d, obj_array);  % Nx x Ny x Nz
	end
	eps_cell = mean_material_node(grid3d, ge, eps_node_cell);
	mu_cell = mean_material_node(grid3d, alter(ge), mu_node_cell);

	% Construct PML s-factors.
	s_factor_cell = generate_s_factor(osc.in_omega0(), grid3d, deg_pml, R_pml);
	pm.mark('eps and mu assignment');

	if ~isTFSF
		% Solve for modes.
		for src = src_array
			if istypesizeof(src, 'ModalSrc')
				modalsrc = src;
				if ~modalsrc.ispreped
					prep_modalsrc(ge, pml, osc, grid3d, eps_cell, mu_cell, s_factor_cell, modalsrc);
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
				cb_center = num2cell(tfsfsrc.shape.cb_center);
				for bgobj = fliplr(obj_array)
					if bgobj.shape.contains(cb_center{:})
						break;  % assume that last object containing TF box center fills TF box
					end
				end
				tfsfsrc.set_bg_material(bgobj.material);
				F0 = tfsfsrc.create_incidentF(osc, grid3d);
				JM = cell(1, Axis.count);
				for w = Axis.elems
					JM{w} = zeros(grid3d.N);
				end
				
				if tfsfsrc.gt == ge
					eqtype_tfsf = EquationType(FT.e, ge);  % for SRCJ, create E-field eq
				else
					eqtype_tfsf = EquationType(FT.h, ge);  % for SRCM, create H-field eq
				end
% 				A = create_eq(eqtype_tfsf, pml, osc.in_omega0(), eps_cell, mu_cell, s_factor_cell, JM, JM, grid3d);
				eq = MatrixEquation(eqtype_tfsf, pml, osc.in_omega0(), eps_cell, mu_cell, s_factor_cell, JM, JM, grid3d);
				Op = eq.matrixfree_op();
							
				x0 = [F0{Axis.x}(:); F0{Axis.y}(:); F0{Axis.z}(:)];
				r = reordering_indices(Axis.count, grid3d.N);
				x0 = x0(r);
							
% 				JM = (A*x0) ./ (-1i*osc.in_omega0());
				JM = Op(x0, 'notransp') ./ (-1i*osc.in_omega0());
				JM = reshape(JM, [Axis.count grid3d.N]);
				JM = permute(JM, [Axis.elems+1, 1]);
				JM = {JM(:,:,:,Axis.x), JM(:,:,:,Axis.y), JM(:,:,:,Axis.z)};
				
				tfsfsrc.setJM(JM, grid3d);
			end
		end
		
		% Add sobj_array to the already-generated eps and mu.
		[eps_node_cell, mu_node_cell] = assign_material_node(grid3d, sobj_array, eps_node_cell, mu_node_cell);  % Nx x Ny x Nz
		eps_cell = mean_material_node(grid3d, ge, eps_node_cell);  % Nx x Ny x Nz
		mu_cell = mean_material_node(grid3d, alter(ge), mu_node_cell);  % Nx x Ny x Nz

		pm.mark('TF/SF source assignment');
	end
	obj_array = [obj_array, sobj_array];
	
	eps_node = cell(1, Axis.count);
	mu_node = cell(1, Axis.count);
	for w = Axis.elems
		eps_node_cell{w} = expand_node_array(grid3d, eps_node_cell{w});  % (Nx+2) x (Ny+2) x (Nz+2)
		mu_node_cell{w} = expand_node_array(grid3d, mu_node_cell{w});  % (Nx+2) x (Ny+2) x (Nz+2)
		
		eps_node{w} = Scalar3d(eps_node_cell{w}, grid3d, [GT.dual GT.dual GT.dual], osc, PhysQ.eps, '\epsilon');
		mu_node{w} = Scalar3d(mu_node_cell{w}, grid3d, [GT.dual GT.dual GT.dual], osc, PhysQ.mu, '\mu');
	end
		
	% Construct sources.
	J_cell = assign_source(grid3d, srcj_array);
	M_cell = assign_source(grid3d, srcm_array);
	pm.mark('J assignment');
end
