%% maxwell_run
% Run MaxwellFDFD.

%%% Syntax
%  [E_cell, H_cell] = maxwell_run(OSC, DOM, OBJ, SRC, [inspect_only])
%  [E_cell, H_cell, obj_array, src_array] = maxwell_run(...)
%  [E_cell, H_cell, obj_array, src_array, extra] = maxwell_run(...)


%%% Description
% |maxwell_run(OSC, DOM, OBJ, SRC)| constructs a simulation domain from given
% objects and sources, and launches a finite-difference frequency-domain (FDFD)
% solver to solve Maxwell's equations in the simulation domain.
%
% Each of |OSC|, |DOM|, |OBJ|, and |SRC| represents a group of parameters.
% Each group supports several flexible expressions.  See the following sections
% about the input parameter groups for more details.
%
% |maxwell_run| can take an optional argument |inspect_only|, which is a logical
% argument (i.e., |true| or |false|).  If |inspect_only = true|, |maxwell_run|
% runs without launching the solver.  This is useful to inspect input arguments
% before starting expensive computation.
% 
% |[E_cell, H_cell] = maxwell_run(...)| returns the E- and H-field solutions of
% Maxwell's equations.  |E_cell| and |H_cell| are length-3 row cell arrays whose
% elements are the x-, y-, and z-components of the respective field solutions.
% Each x-, y-, z-component of the fields are provided as instances of
% <Scalar3d.html Scalar3d>.
%
% |[E_cell, H_cell, obj_array, src_array] = maxwell_run(...)| returns arrays of
% instances of <Object.html |Object|> and <Source.html |Source|>.  The |Object|
% and |Source| elements represent the objects and sources placed in the
% simulation domain, so they can be used to visualize the simulation domain.
%
% |[E_cell, H_cell, obj_array, src_array, extra] = maxwell_run(...)| returns
% |extra|, which is a structure containing the fields |grid3d|, |J|, |M|, |eps|,
% and |mu|: |grid3d| is an instance of <Grid3d.html Grid3d> that contains the
% grid information; |J| and |M| are the electric and magnetic current source
% distributions that have the same structure as |E_cell| and |H_cell|; |eps| and
% |mu| are instances of <Scalar3d.html Scalar3d> that contain the electric
% permittivity and magnetic permeability defined at the centers of the grid
% cells.


%%% Input Parameter Group - OSC
% OSC group specifies oscillation parameters.  The group begins with |'OSC'| and
% ends with one of the followings:
%
% * |L0, wvlen|  
% * |osc|
%
% |L0|: unit of length in meters. For example, |L0 = 1e-9| means 1 nm.
% All lengths in other parameters are described in this unit. 
%
% |wvlen|: vacuum wavelength in the unit of |L0|.  In other words, the actual
% wavelength is |wvlen * L0|.
%
% |osc|: instance of <Oscillation |Oscillation|>, which contains the information
% about the length unit and vacuum wavelength.


%%% Input Parameter Group - DOM
% DOM group specifies domain parameters.  The group begins with |'DOM'| and ends
% with one of the followings:
%
% * |material, box, dl, bc, Lpml, [deg_pml, [R_pml]], [withuniformgrid]|
% * |domain, dl, bc, Lpml, [deg_pml, [R_pml]], [withuniformgrid]|
%
% |material|: background material filling the simulation domain.  See
% <maxwell_run.html#7 Material Description> for various ways to describe a
% material.
%
% |box|: range of the simulation doamin in the form of |[xmin xmax; ymin ymax;
% zmin zmax]|.
%
% |dl|: default grid cell size in the form of |[dx dy dz]|.  Can be abbreviated
% to a scalar |dl| if |dx == dy == dz|.
%
% |bc|: boundary conditions in the form of |[bc_xn bc_yn bc_zn]|, whose each
% element specifies the boundary condition at the negative end in one of the x-,
% y-, z-axes (e.g., |bc_xn| for the negative end in the x-axis.).  Each element
% of |bc| is an instance of <BC.html |BC|>.  The boundary conditions at the
% positive ends are automatically determined by those at the negative ends:
% |bc_wp == BC.p| if |bc_wn == BC.p|, and |bc_wp == BC.m| otherwise.  Can be
% further abbreviated to a single BC instance if |bc_xn == bc_yn == bc_zn|.
%
% |Lpml|: thicknesses of PML in the form of |[Lpml_xn Lpml_xp; Lpml_yn Lpml_yp;
% Lpml_zn Lpml_zp]|, whose each element specifies the thickness at either
% negative or positive end in one of the x-, y-, z-axes (e.g., |Lpml_xn| for the
% negative end in the x-axis).  Can be abbreviated to |[Lpml_x Lpml_y Lpml_z]|
% if |Lpml_wn == Lpml_wp| for |w = x, y, z|.  Can be further abbreviated to a
% single scalar thickness if |Lpml_x == Lpml_y == Lpml_z|.  All thicknesses are
% in the unit of |L0|.
%
% |domain|: instance of <Domain.html |Domain|>, which contains the information
% about |material| and |box|.
%
% |deg_pml|: optional argument to specify the degrees of the polynomial gradings
% of PML loss parameters in the form of |[deg_xn deg_xp; deg_yn deg_yp; deg_zn
% deg_zp]|.  Can be abbreviated to |[deg_x deg_y deg_z]| if |deg_wn == deg_wp|
% for |w = x, y, z|.  Can be further abbreviated to a single scalar if |deg_x ==
% deg_y == deg_z|.  The default value is |deg_pml = 4|, meaning that the
% polynomials of degree 4 are used to smoothly increase the PML loss parameters
% for all boundaries with PML.
%
% |R_pml|: optional argument to specify the target reflectances in the form of
% |[R_xn R_xp; R_yn R_yp; R_zn R_zp]|.  To specify |R_pml|, |deg_pml| should be
% specified.  Can be abbreviated to |[R_x R_y R_z]| if |R_wn == R_wp| for |w =
% x, y, z|.  Can be further abbreviated to a single scalar if |R_x == R_y ==
% R_z|.  The default value is |R_pml = exp(-16)|, meaning that the reflectance
% from PML is aimed to be as small as |exp(-16) ~= 1e-7|.
%
% |withuniformgrid|: optional argument to select a grid generation algorithm.
% If |withuniformgrid = true|, a grid is generated uniformly; otherwise a grid
% is generated dynamically to best fit the objects in the simulation domain.
% The default value is |withuniformgrid = false|, i.e., dynamic grid generation
% is preferred to uniform grid generation.


%%% Input Parameter Group - OBJ and SOBJ
% OBJ group specifies obect parameters.  The group begins with |'OBJ'| and ends
% with one of the followings:
%
% * |material_1, shapes_1, ..., material_N, shapes_N|
% * |obj_1, ..., obj_N|
% * |eps_array|
%
% |material_i|: material filling |shapes_i|.  See <maxwell_run.html#7 Material
% Description> for various ways to describe a material.
%
% |shapes_i|: shapes made of |material_i|.  It can be a single shape, an array
% of shapes (i.e., |[shape_a, shape_b, ...]|), or can be spanned into several
% shape arguments (i.e., |shape_1, shape_2, ...|).  Each shape is an instance of
% <Shape.html |Shape|>.
%
% |obj_i|: instance or array of instances of <Object.html |Object|>.  Each
% |Object| is composed of a material and shape.
%
% |eps_array| is a 3D array of permittivities defined at grid vertices.  The
% size of the array should be the same as the number of vertices in a generated
% grid.  Because the number of vertices is hard to calculate for dynamic grid
% generation, uniform grid generation is recommended when |eps_node_array| is
% used.
%
% There is a similar, optional parameter group SOBJ that begins with |'SOBJ'|.
% SOBJ stands for _scattering objects_, and it is used to define scatterers for
% total-field/scattered-field (TF/SF) simulation.  When a TF/SF source is used
% as a source, the objects defined in SOBJ group are treated as scatterers
% whereas the objects defined in OBJ group are treated as background objects.
% If a TF/SF source is not used, SOBJ group works the same as OBJ group.


%%% Input Parameter Group - SRC
% SRC group specifies current source parameters. Depending on whether the
% sources are electric or magnetic, the group begins with |'SRCJ'| or |'SRCM'|
% and ends with
%
% * |src_1, ..., src_M|
%
% |src_i|: instance or array of instances of <Source.html |Source|>.  


%%% Material Description
% Each material is described by one of the followings:
%
% * |{name, color, permittivity}|
% * |{material_datapath, color}|
% * |material|
%
% |name|: string describing the name of the material (e.g., 'vacuum', 'silver',
% 'Ag'). 
%
% |color|: color to be used in visualizing objects made of the material.  It can
% be either a standard color specifier character used in MATLAB (e.g., |'w'| for
% white, |'k'| for black, etc.) or |[r, g, b]| where |r|, |g|, |b| are RGB color
% scales ranging from |0.0| to |1.0| (e.g., |[0.5, 0.5, 0.5]| for gray).  Use
% |color = 'none'| to not draw the material.
%
% |permittivity|: electric permittivity value in the unit of the vacuum
% permittivity.  It is a complex number.
%
% |material_datapath|: path to the file of the permittivity data of the
% material.  It is a path relative to |material| directory.  For example,
% |'Palik/Si'| refers to silicon's permittivity data stored in |Si.m| file in
% |material/Palik| directory.  The permittivity for the frequency of interest is
% automatically retrieved from the data file.  Note that |"/"| here is the file
% separator used in UNIX-like systems; In MS Windows, |"\"| should be used.  To
% be platform-independent, |filesep| can be used in place of |"/"| or |"\"|
% (e.g., |strcat('Palik', filesep, 'Si')| or |['Palik', filesep, 'Si']|).


%%% Example
%   gray = [0.5 0.5 0.5];  % [r g b]
%   inspect_only = true;
%   [E, H, obj_array, extra] = maxwell_run(...
%       'OSC', 1e-9, 1550, ...
%       'DOM', {['Palik', filesep, 'SiO2'], 'none'}, [-700, 700; -600, 600; -200, 1700], 20, BC.p, 200, ...
%       'OBJ', ...
%           {['Palik', filesep, 'SiO2'], 'none'}, Box([-50, 50; -50, 50; -200, 1700], [2, 2, 20]), ...  % OBJ1
%           {['CRC', filesep, 'Ag'], gray}, [Box([-700, -25; -25, 25; -200, 1700], 20), Box([25, 700; -25, 25; -200, 1700], 20)], ...  % OBJ2
%       'SRCJ', PointSrc(Axis.x, [0, 0, 200]), ...
%       inspect_only);

function [E_cell, H_cell, obj_array, src_array, extra] = maxwell_run(varargin)
	DEFAULT_METHOD = 'direct';  % 'direct', 'gpu', 'aws', 'inputfile'
		
	% Set solver options.
	iarg = nargin; arg = varargin{iarg};
	inspect_only = false;
	if istypesizeof(arg, 'logical')
		inspect_only = arg;
		iarg = iarg - 1; arg = varargin{iarg};
	end

	is_solveropts = false;
	if istypesizeof(arg, 'struct')
		solveropts = arg;
		is_solveropts = true;
		iarg = iarg - 1;  % arg = varargin{iarg};
	end
		
	if ~is_solveropts || ~isfield(solveropts, 'showstruct')
		solveropts.showstruct = true;
	end
	
	if ~is_solveropts || ~isfield(solveropts, 'method')
		solveropts.method = DEFAULT_METHOD;
	end
	
	if is_solveropts && isequal(solveropts.method, 'inputfile')
		chkarg(isfield(solveropts, 'filenamebase'), '"solveropts" should have "filenamebase" field.');
	end
	
	if ~is_solveropts || ~isfield(solveropts, 'maxit')
% 		solveropts.maxit = intmax;
		solveropts.maxit = 1e6;
	else
		chkarg(istypesizeof(solveropts.maxit, 'real') && solveropts.maxit > 0, ...
			'solveropts.maxit should be positive.');	
	end

	if ~is_solveropts || ~isfield(solveropts, 'tol')
		solveropts.tol = 1e-6;
	else
		chkarg(istypesizeof(solveropts.tol, 'real') && solveropts.tol > 0, ...
			'solveropts.tol should be positive.');
	end

	if ~is_solveropts || ~isfield(solveropts, 'eqtype')
		solveropts.eqtype = EquationType(FT.e, GT.prim);
	else
		chkarg(istypesizeof(solveropts.eqtype, 'EquationType'), ...
			'solveropts.eqtype should be instance of EquationType.');
	end

	if ~is_solveropts || ~isfield(solveropts, 'pml')
		solveropts.pml = PML.sc;
	else
		chkarg(istypesizeof(solveropts.pml, 'PML'), ...
			'solveropts.pml should be instance of PML.');
	end
	
	if ~is_solveropts || ~isfield(solveropts, 'returnAandb')
		solveropts.returnAandb = false;
	else
		chkarg(istypesizeof(solveropts.returnAandb, 'logical'), ...
			'solveropts.returnAandb should be logical.');
	end
	
	if ~is_solveropts || ~isfield(solveropts, 'returnDiv')
		solveropts.returnDiv = false;
	else
		chkarg(istypesizeof(solveropts.returnDiv, 'logical'), ...
			'solveropts.returnDiv should be logical.');
	end
	
	chkarg(iarg > 0, 'first argument is not correct.');

	if inspect_only
		fprintf('%s begins (inspection only).\n', mfilename);
	else
		fprintf('%s begins.\n', mfilename);
	end
	ge = solveropts.eqtype.ge;
	fprintf('E-field grid type: %s\n', char(ge));
	pm = ProgMark();
	
	% Build the system.
	% Make sure to pass the first consecutive elements of varargin to
	% build_system() for correct error reports.
	[osc, grid3d, s_factor, eps, mu, J, M, obj_array, src_array, mat_array, eps_node, mu_node, isiso] = ...
		build_system(ge, solveropts.pml, varargin{1:iarg}, pm);
	
	if ~is_solveropts || ~isfield(solveropts, 'F0')
		solveropts.F0 = 'zero';  % 'rand' is the other choice
	else
		chkarg(isequal(solveropts.F0, 'zero') || isequal(solveropts.F0, 'rand') ...
			|| istypesizeof(solveropts.F0, 'complexcell', [1 Axis.count], grid3d.N), ...
			'solveropts.F0 should be length-%d cell array whose each element is %d-by-%d-by-%d array with complex numbers.', ...
			Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
	end
		
	if inspect_only  % inspect objects and sources
		if solveropts.showstruct
			figure;
			set(gcf, 'units','normalized','position',[0.5 0 0.5 0.5]);			
			withpml = true;
			visobjsrc(grid3d, obj_array, src_array, withpml);
			drawnow
			pm.mark('domain visualization');
		end
		
		% Visualize modes.
		is_modalsrc = false;
		for src = src_array
			if istypesizeof(src, 'ModalSrc')
				is_modalsrc = true;
				modalsrc = src;
				opts.withabs = true;
				
				E2d = modalsrc.E2d;
				cmax = max(abs([E2d{Axis.x}.array(:); E2d{Axis.y}.array(:); E2d{Axis.z}.array(:)]));
				opts.cmax = cmax;
				for w = Axis.elems
					figure;
					set(gcf, 'units','normalized','position',[(int(w)-1)/3 1/2 1/3 1/3]);			
					vis2d(E2d{w}, obj_array, opts);
					drawnow;
				end
				
				H2d = modalsrc.H2d;
				cmax = max(abs([H2d{Axis.x}.array(:); H2d{Axis.y}.array(:); H2d{Axis.z}.array(:)]));
				opts.cmax = cmax;
				for w = Axis.elems
					figure;
					set(gcf, 'units','normalized','position',[(int(w)-1)/3 0 1/3 1/3]);			
					vis2d(H2d{w}, obj_array, opts);
					drawnow;
				end
			end
		end
		
		if is_modalsrc
			pm.mark('distributed source visualization');
		end
		fprintf('%s finishes (inspection only).\n\n', mfilename);
		E = {};
		H = {};
	else  % inspect_only == false
		if isequal(solveropts.method, 'inputfile')
			% Define eps_node_array at vertices of the E-field edges.
			chkarg(isiso, 'anisotropic materials are not supported for solveropts.method == ''inputfile''.');
			
			if ge == GT.prim
				eps_node = eps_node{Axis.x};  % ad hoc solution
				
				% Old implementation (without anisotropy support) starts
				% here.
				eps_node_array = eps_node.data_expanded();  % (Nx+2) x (Ny+2) x (Nz+2)
				eps_node_array = eps_node_array(1:end-1, 1:end-1, 1:end-1);  % (Nx+1) x (Ny+1) x (Nz+1)
				
				% The below line is an ad hoc solution for the error.  It
				% is just to pass an array with correct size to
				% write_input(), but eps_node_array is not used in
				% write_input() currently.
				eps_node_array = eps_node_array(1:end-1, 1:end-1, 1:end-1);  % Nx x Ny x Nz
% 
% 				% (Anisotropy support? begins)
% % 				eps_node_cell = eps_cell;
% % 				for w = Axis.elems
% % 					ind_g = {':',':',':'};
% % 					if grid3d.bc(w) == BC.p
% % 						ind_g{w} = grid3d.N(w);
% % 					else
% % 						ind_g{w} = 1;
% % 					end
% % 					eps_node_cell{w} = cat(int(w), eps_node_cell{w}(ind_g{:}), eps_node_cell{w});
% % 				end
% % 				eps_node_cell{Axis.x} = 2./(1./eps_node_cell{Axis.x}(1:end-1, :, :) + 1./eps_node_cell{Axis.x}(2:end, :, :));
% % 				eps_node_cell{Axis.y} = 2./(1./eps_node_cell{Axis.y}(:, 1:end-1, :) + 1./eps_node_cell{Axis.y}(:, 2:end, :));
% % 				eps_node_cell{Axis.z} = 2./(1./eps_node_cell{Axis.z}(:, :, 1:end-1) + 1./eps_node_cell{Axis.z}(:, :, 2:end));
% % 				
% % 				eps_node_array = (eps_node_cell{Axis.x} + eps_node_cell{Axis.y} + eps_node_cell{Axis.z})./3;
% % 	
% % 				eps_node_array = (eps_node_array(1:end-1,1:end-1,1:end-1) + eps_node_array(2:end,1:end-1,1:end-1) ...
% % 							+ eps_node_array(1:end-1,2:end,1:end-1) + eps_node_array(1:end-1,1:end-1,2:end) ...
% % 							+ eps_node_array(1:end-1,2:end,2:end) + eps_node_array(2:end,1:end-1,2:end) ...
% % 							+ eps_node_array(2:end,2:end,1:end-1) + eps_node_array(2:end,2:end,2:end))./8;
% 				% (Anisotropy support? ends)
% 
% 				eps_node_array = 8./(1./eps_node_array(1:end-1,1:end-1,1:end-1) + 1./eps_node_array(2:end,1:end-1,1:end-1) ...
% 							+ 1./eps_node_array(1:end-1,2:end,1:end-1) + 1./eps_node_array(1:end-1,1:end-1,2:end) ...
% 							+ 1./eps_node_array(1:end-1,2:end,2:end) + 1./eps_node_array(2:end,1:end-1,2:end) ...
% 							+ 1./eps_node_array(2:end,2:end,1:end-1) + 1./eps_node_array(2:end,2:end,2:end));
			else
				eps_node_array = eps_node.data_original();  % Nx x Ny x Nz
			end

			write_input(solveropts.filenamebase, solveropts.eqtype, osc, grid3d, s_factor, ...
				eps_node_array, eps, mu, J, M, solveropts.F0, solveropts.tol, solveropts.maxit);

			pm.mark('input file creation');		
			fprintf('%s finishes. (input file created)\n\n', mfilename);
			E = {};
			H = {};
		else  % solveropts.method ~= 'inputfile'
			if isequal(solveropts.method, 'direct')
				[E, H] = solve_eq_direct(solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
			elseif isequal(solveropts.method, 'iterative')
				[E, H, relres, iter, resvec] = solve_eq_iterative(solveropts.maxit, solveropts.tol, solveropts.F0, solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
			else  % for solvers using E-field based equation
				%% Apply spatial inversion.
% 				d_prim = grid3d.dl(:, GT.prim);
% 				d_dual = grid3d.dl(:, GT.dual);
% 				s_prim = s_factor(:, GT.prim);
% 				s_dual = s_factor(:, GT.dual);
				d_prim = flip_vec(grid3d.dl(:, GT.dual));  % GT.dual, not GT.prim
				d_dual = flip_vec(grid3d.dl(:, GT.prim));  % GT.prim, not GT.dual
				s_prim = flip_vec(s_factor(:, GT.dual));  % GT.dual, not GT.prim
				s_dual = flip_vec(s_factor(:, GT.prim));  % GT.prim, not GT.dual
				mu = flip_vec(mu);
				eps = flip_vec(eps);
				J = neg_vec(flip_vec(J));  % pseudovector

				if isequal(solveropts.method, 'gpu')
					ds_prim = mult_vec(d_prim, s_prim);
					ds_dual = mult_vec(d_dual, s_dual);
					figure;
					F0 = {zeros(grid3d.N), zeros(grid3d.N), zeros(grid3d.N)};
					[E, H] = fds(osc.in_omega0(), ...
									ds_prim, ds_dual, ...
									mu, eps, ...
									F0, J, ...
									solveropts.maxit, solveropts.tol, 'plot');
					%   norm(A2 * ((1./e) .* (A1 * y)) - omega^2 * m .* y - A2 * (b ./ (-i*omega*e))) / norm(b) % Error for H-field wave equation.
				elseif isequal(solveropts.method, 'aws')
					ds_prim = mult_vec(d_prim, s_prim);
					ds_dual = mult_vec(d_dual, s_dual);
					F0 = {zeros(grid3d.N), zeros(grid3d.N), zeros(grid3d.N)};
% 					F0 = {rand(grid3d.N), rand(grid3d.N), rand(grid3d.N)};
% 					F0 = {rand(1)*ones(grid3d.N), rand(1)*ones(grid3d.N), rand(1)*ones(grid3d.N)};
					callback = maxwell(osc.in_omega0(), ...
									ds_prim, ds_dual, ...
									mu, eps, ...
									F0, J, ...
									solveropts.maxit, solveropts.tol);
					while ~callback(); end
					[~, E, H] = callback();
				end

				E = neg_vec(flip_vec(E));  % pseudovector
				J = neg_vec(flip_vec(J));  % pseudovector
				H = flip_vec(H);
			end

			pm.mark('solution calculation');
			if isequal(solveropts.method, 'iterative')
				%% Report the iterative solution process.
				semilogy(1:length(resvec), resvec);
				fprintf('\t# of iteration steps = %d\n', iter);
				fprintf('\trelative residual error = %e\n', relres);
			end
			fprintf('unknowns: %s-field\n', char(solveropts.eqtype.f));
			fprintf('%s finishes.\n\n', mfilename);
		end
	end
	
	if solveropts.returnAandb
%		[A, b, g_from_f] = create_eq(solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
		eq = MatrixEquation(solveropts.eqtype, solveropts.pml, osc.in_omega0(), eps, mu, s_factor, J, M, grid3d);
		
		[extra.A, extra.b] = eq.matrix_op();
		[~, ~, extra.GfromF] = eq.matrixfree_op();
	end

	if solveropts.returnDiv
		[Dive, Divm] = create_divs(ge, s_factor, grid3d);
		extra.Dive = Dive;
		extra.Divm = Divm;
	end

	% Construct Scalar3d objects.
	E_cell = cell(1, Axis.count);
	H_cell = cell(1, Axis.count);
	J_cell = cell(1, Axis.count);
	M_cell = cell(1, Axis.count);
	for w = Axis.elems
		if ~isempty(E)
			E_cell{w} = array2scalar(E{w}, grid3d, ge, FT.e, w, osc, PhysQ.E);
		end
		
		if ~isempty(H)
			H_cell{w} = array2scalar(H{w}, grid3d, ge, FT.h, w, osc, PhysQ.H);
		end
		
		J_cell{w} = array2scalar(J{w}, grid3d, ge, FT.e, w, osc, PhysQ.J);
		M_cell{w} = array2scalar(M{w}, grid3d, ge, FT.h, w, osc, PhysQ.M);
	end
	
	extra.grid3d = grid3d;
	extra.J = J_cell;
	extra.M = M_cell;
	extra.eps = eps_node;
	extra.mu = mu_node;
end
