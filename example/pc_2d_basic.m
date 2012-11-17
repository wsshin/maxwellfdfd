clear all; close all; clear classes; clc;

%% Set flags.
isnew = true;
inspect_only = false;

%%
if isnew
	%% Create shapes.
	a = 420;  % lattice constant
	t = 1;  % slab thickness
	r = 0.29*a;  % hole radius

	ad = 25;  % divider for a
	td = 10;  % divider for t
	dd = 10;  % divider for d = 2*r

	mx = 7.5;
	my = 4.5;
	slab_yn = Box([-mx*a mx*a; -my*a -0.5*a; 0 t], [a/ad, a/ad, t]);
	slab_yp = Box([-mx*a mx*a; 0.5*a my*a; 0 t], [a/ad, a/ad, t]);

	rod = CircularCylinder(Axis.z, [0 0 t/2], r, t, [2*r/dd, 2*r/dd, t]);

	%% Solve the system.
	gray = [0.5 0.5 0.5];  % [r g  b]
	solveropts.method = 'direct';

% 	[E, H, obj_array, err] = maxwell_run(c, 1e-9, 1550, ...
	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1550, ...
		{'vacuum', 'white', 1.0}, [-mx*a mx*a; -my*a my*a; 0 t], [a/ad a/ad t], BC.p, [2*a 2*a 0], ...
		{'dieletric', gray, 11.56}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yn), ...
		{'dielectric', gray, 11.56}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yp), ...
		PointSrc(Axis.z, [0, 0, 0]), solveropts, inspect_only);

	if ~inspect_only
		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end


%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts
% 	opts.cscale = 5e-3;
% 	opts.withobj = true;
	opts.withobj = true;
%  	opts.withgrid = true;
% 	opts.withinterp = false;
% 	visall(E{Axis.y}, obj_array, opts);
	vis2d(E{Axis.z}, Axis.z, 0, obj_array, src_array, opts)  % why don't we duplicate ezplot.m and modify it so that it uses a finer mesh?

% %%
% 	Sx = poynting(Axis.x, E{Axis.y}, E{Axis.z}, H{Axis.y}, H{Axis.z}, Axis.x, 100);
% 	figure;
% 	vis2d(Sx, obj_array);
% 
% 	p = flux_patch(Sx)
end
