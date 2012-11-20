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

	mx = 20.5;
	my = 7.5;
	slab_yn = Box([-mx*a mx*a; -my*a -0.5*a; 0 t], [a/ad, a/ad, t]);
	slab_yp = Box([-mx*a mx*a; 0.5*a my*a; 0 t], [a/ad, a/ad, t]);

	rod = CircularCylinder(Axis.z, [0 0 t/2], r, t, [2*r/dd, 2*r/dd, t]);

	%% Solve the system.
	gray = [0.5 0.5 0.5];  % [r g  b]
	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1550, ...
		{'vacuum', 'white', 1.0}, [-mx*a mx*a; -my*a my*a; 0 t], [a/ad a/ad t], BC.p, [5*a 2*a 0], ...
		{'dieletric', gray, 11.56}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yn), ...
		{'dielectric', gray, 11.56}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yp), ...
		PointSrc(Axis.z, [0, 0, 0]), inspect_only);

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
% 	opts.withobj = true;
	opts.withabs = true;
	opts.withpml = false;
	vis2d(E{Axis.z}, Axis.z, 0, obj_array, src_array, opts)
end
