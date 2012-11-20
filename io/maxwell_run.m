%% maxwell_run
% Run MaxwellFDS.

%%% Syntax
%  [E_cell, H_cell] = maxwell_run(FREQ, DOMAIN, OBJ, SRC, [inspect_only])
%  [E_cell, H_cell, obj_array, src_array] = maxwell_run(...)
%  [E_cell, H_cell, obj_array, src_array, err] = maxwell_run(...)


%%% Description
% |maxwell_run(FREQ, DOMAIN, OBJ, SRC)| constructs a simulation domain from
% given objects and sources, and launches the frequency-domain solver (FDS) to
% solve Maxwell's equations in the simulation domain.
%
% |maxwell_run| can take an optional argument |inspect_only|, which is a logical
% argument (i.e., |true| or |false|).  If |inspect_only = true|, |maxwell_run|
% runs without launching FDS.  This is useful to inspect input arguments before
% starting expensive computation.
% 
% Each of |FREQ|, |DOMAIN|, |OBJ|, and |SRC| represents a group of parameters.
% Each group supports several flexible expressions.  See the following sections
% about the input parameter groups for more details.
%
% |[E_cell, H_cell] = maxwell_run(...)| returns the E- and H-field solutions of
% Maxwell's equations.  |E_cell| and |H_cell| are length-3 row cell arrays whose
% elements are the x-, y-, and z-components of the respective field solutions.
% Each x-, y-, z-component of the fields are provided as instances of
% <Scalar3d.html Scalar3d>.

% |[E_cell, H_cell, obj_array, src_array] = maxwell_run(...)| returns arrays of
% instances of <Object.html |Object|> and <Source.html |Source|>.  The |Object|
% and |Source| elements represent the objects and sources placed in the
% simulation domain, so they can be used to visualize the simulation domain.

% |[E_cell, H_cell, obj_array, src_array, err] = maxwell_run(...)| returns a
% column vector |err| that records the relative residual errors in the field
% solutions as the solutions are evolved by the iterative method used in FDS.


%%% Input Parameter Group - FREQ
% Possible expressions of FREQ group are listed below:
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


%%% Input Parameter Group - DOMAIN
% Possible expressions of DOMAIN group are listed below:
%
% * |material, box, dl, bc, Lpml, [withuniformgrid]|
% * |domain, dl, bc, Lpml, [withuniformgrid]|
%
% |material|: background material filling the simulation domain.  See
% <maxwell_run.html#7 Material Description> for various ways to describe a
% material.
%
% |box|: range of the simulation doamin in the form of |[xmin xmax; ymin ymax;
% zmin zmax]|.
%
% |dl|: default grid cell size in the form of |[dx, dy, dz]|.  Can be
% abbreviated to a scalar |dl| if |dx == dy == dz|.
%
% |bc|: boundary conditions in the form of |[bc_xn, bc_yn, bc_zn]|, whose each
% element specifies the boundary condition at the negative end in one of the x-,
% y-, z-axes (e.g., |bc_xn| for the negative end in the x-axis.).  Each element
% of |bc| is an instance of <BC.html |BC|>.  The boundary conditions at the
% positive ends automatically determines by those at the negative ends: |bc_wp
% == BC.p| for |bc_wn == BC.p| and |bc_wp == BC.m| otherwise.  Can be further
% abbreviated to a single BC instance if |bc_xn == bc_yn == bc_zn|.
%
% |Lpml|: thicknesses of PML in the form of |[Lpml_xn, Lpml_xp; Lpml_yn,
% Lpml_yp; Lpml_zn, Lpml_zp]|, whose each element specifies the thickness at
% either negative or positive end in one of the x-, y-, z-axes (e.g., |Lpml_xn|
% for the negative end in the x-axis).  Can be abbreviated to |[Lpml_x, Lpml_y,
% Lpml_z]| if |Lpml_wn == Lpml_wp| for |w = x, y, z|.  Can be further
% abbreviated to a single scalar thickness if |Lpml_x == Lpml_y == Lpml_z|.  All
% thicknesses are in the unit of |L0|.
%
% |domain|: instance of <Domain.html |Domain|>, which contains the information
% about |material| and |box|.
%
% |withuniformgrid|: optional argument to select a grid generation algorithm.
% If |withuniformgrid = true|, a grid is generated uniformly; otherwise a grid
% is generated dynamically to best fit the objects in the simulation domain.
% The default value is |withuniformgrid = false|, i.e., dynamic grid generation
% is preferred to uniform grid generation.


%%% Input Parameter Group - OBJ
% Possible expressions of OBJ group are listed below:
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


%%% Input Parameter Group - SRC
% Possible expressions of DOMAIN group are listed below:
%
% * |src_1, ..., src_M|
%
% |src_i|: instance or array of instances of <Source.html |Source|>.  


%%% Material Description
% Possible expressions for each material description are listed below:
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
%   [E, H, obj_array, err] = maxwell_run(1e-9, 1550, ...  % FREQ
%       {['Palik', filesep, 'SiO2'], 'none'}, [-700, 700; -600, 600; -200, 1700], 20, BC.p, 200, ...  % DOMAIN
%       {['Palik', filesep, 'SiO2'], 'none'}, Box([-50, 50; -50, 50; -200, 1700], [2, 2, 20]), ...  % OBJ1
%       {['CRC', filesep, 'Ag'], gray}, [Box([-700, -25; -25, 25; -200, 1700], 20), Box([25, 700; -25, 25; -200, 1700], 20)], ...  % OBJ2
%       PointSrc(Axis.x, [0, 0, 200]), ...  % SRC
%       inspect_only ...
%       );

function [E_cell, H_cell, obj_array, src_array, err] = maxwell_run(varargin)
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
		
	if ~is_solveropts || ~isfield(solveropts, 'method')
		solveropts.method = DEFAULT_METHOD;
	end
	
	if is_solveropts && isequal(solveropts.method, 'aws')
		chkarg(isfield(solveropts, 'cluster') && isfield(solveropts, 'nodes'), ...
			'"solveropts" should have "cluster" and "nodes" fields.');
	end

	if is_solveropts && isequal(solveropts.method, 'inputfile')
		chkarg(isfield(solveropts, 'filenamebase'), '"solveropts" should have "filenamebase" field.');
	end
	
	if ~is_solveropts || ~isfield(solveropts, 'maxit')
		solveropts.maxit = 1e5;
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

	chkarg(iarg > 0, 'first argument is not correct.');

	if inspect_only
		fprintf('%s begins (inspection only).\n', mfilename);
	else
		fprintf('%s begins.\n', mfilename);
	end
	pm = ProgMark();
	
	% Build the system.
	% Make sure to pass the first consecutive elements of varargin to
	% build_system() for correct error reports.
	[osc, grid3d, s_factor, eps_face, mu_edge, J, E0, obj_array, src_array, eps_node] = ...
		build_system(varargin{1:iarg}, pm);
	
	if inspect_only  % inspect objects and sources
		figure;
		set(gcf, 'units','normalized','position',[0.5 0 0.5 0.5]);			
		withpml = false;
		visobjsrc(grid3d, obj_array, src_array, withpml);
		drawnow
		pm.mark('domain visualization');
		
		% Solve for modes.
		is_distsrc = false;
		for src = src_array
			if istypesizeof(src, 'DistributedSrc')
				is_distsrc = true;
				distsrc = src;
				[h, v, n] = cycle(distsrc.normal_axis);
				grid2d = Grid2d(grid3d, n);

				Jh2d = array2scalar(distsrc.Jh, PhysQ.J, grid2d, h, GK.dual, osc, distsrc.intercept);
				Jv2d = array2scalar(distsrc.Jv, PhysQ.J, grid2d, v, GK.dual, osc, distsrc.intercept);

				cmax = max(abs([Jh2d.array(:); Jv2d.array(:)]));
				opts.withabs = true;
				opts.cmax = cmax;
				figure;
				set(gcf, 'units','normalized','position',[0 0.5 0.5 0.5]);			
				vis2d(Jh2d, obj_array, opts);
				drawnow;

				figure;
				set(gcf, 'units','normalized','position',[0.5 0.5 0.5 0.5]);			
				vis2d(Jv2d, obj_array, opts);
				drawnow;
			end
		end
		
		if is_distsrc
			pm.mark('distributed source visualization');
		end
	
		E_cell = {};
		H_cell = {};
		fprintf('%s finishes (inspection only).\n\n', mfilename);
		return;
	elseif isequal(solveropts.method, 'inputfile')
		write_input(solveropts.filenamebase, osc, grid3d, s_factor, ...
			eps_node(1:end-1,1:end-1,1:end-1), eps_face, mu_edge, J, solveropts.tol, solveropts.maxit);

		E_cell = {};
		H_cell = {};
		pm.mark('input file creation');		
		fprintf('%s finishes. (input file created)\n\n', mfilename);
		return;
	end
	
	%% Apply spatial inversion.
    flip_array = @(F) flipdim(flipdim(flipdim(F, int(Axis.z)), int(Axis.y)), int(Axis.x));
    flip_vec = @(F_cell) {flip_array(F_cell{Axis.x}), flip_array(F_cell{Axis.y}), flip_array(F_cell{Axis.z})};
    neg_vec = @(F_cell) {-F_cell{Axis.x}, -F_cell{Axis.y}, -F_cell{Axis.z}};
	mult_vec = @(F_cell, G_cell) {F_cell{Axis.x}.*G_cell{Axis.x}, F_cell{Axis.y}.*G_cell{Axis.y}, F_cell{Axis.z}.*G_cell{Axis.z}};
		
% 	d_prim = grid3d.dl(:, GK.prim);
% 	d_dual = grid3d.dl(:, GK.dual);
% 	s_prim = s_factor(:, GK.prim);
% 	s_dual = s_factor(:, GK.dual);
	d_prim = flip_vec(grid3d.dl(:, GK.dual));  % GK.dual, not GK.prim
	d_dual = flip_vec(grid3d.dl(:, GK.prim));  % GK.prim, not GK.dual
	s_prim = flip_vec(s_factor(:, GK.dual));  % GK.dual, not GK.prim
	s_dual = flip_vec(s_factor(:, GK.prim));  % GK.prim, not GK.dual
	mu_edge = flip_vec(mu_edge);
	eps_face = flip_vec(eps_face);
	J = neg_vec(flip_vec(J));  % pseudovector
	E0 = neg_vec(flip_vec(E0));  % pseudovector
	
	if isequal(solveropts.method, 'direct')
		[E, H] = solve_eq_direct(osc.in_omega0(), ...
						d_prim, d_dual, ...
						s_prim, s_dual, ...
						mu_edge, eps_face, ...
						J, E0);
	elseif isequal(solveropts.method, 'gpu')
		ds_prim = mult_vec(d_prim, s_prim);
		ds_dual = mult_vec(d_dual, s_dual);
		figure;
%		[E, H, err] = solve(cluster, osc.in_omega0(), ...
		[E, H, err] = fds(osc.in_omega0(), ...
						ds_prim, ds_dual, ...
						mu_edge, eps_face, ...
						E0, J, ...
						solveropts.maxit, solveropts.tol, 'plot');
		%   norm(A2 * ((1./e) .* (A1 * y)) - omega^2 * m .* y - A2 * (b ./ (-i*omega*e))) / norm(b) % Error for H-field wave equation.
	elseif isequal(solveropts.method, 'aws')
		[E, H, err] = maxwell.solve(solveropts.cluster, solveropts.nodes, ...
						osc.in_omega0(), ...
						d_prim, d_dual, ...
						s_prim, s_dual, ...
						mu_edge, eps_face, ...
						E0, J, ...
						solveropts.maxit, solveropts.tol, 'plot');
	end
	
	E = neg_vec(flip_vec(E));  % pseudovector
	H = flip_vec(H);
	
	pm.mark('solution calculation');

	% Construct Scalar3d objects.
	E_cell = cell(1, Axis.count);
	H_cell = cell(1, Axis.count);
	for w = Axis.elems
		E_cell{w} = array2scalar(E{w}, PhysQ.E, grid3d, w, GK.dual, osc);
		H_cell{w} = array2scalar(H{w}, PhysQ.H, grid3d, w, GK.prim, osc);
	end
	pm.mark('solution preparation');
		
	fprintf('%s finishes.\n\n', mfilename);
end
