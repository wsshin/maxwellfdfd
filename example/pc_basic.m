clear all; close all; clear classes; clc;

%%
isnew = false;
% inspect_only = true;
inspect_only = false;

%% Create shapes.
if isnew
	a = 420;  % lattice constant
	t = 0.6*a;  % slab thickness
	r = 0.29*a;  % hole radius
	h = sqrt(3)/2*a;  % distance between rows of holes

	ad = 20;  % divider for a
	td = 10;  % divider for t
	dd = 10;  % divider for d = 2*r

	slab = Box([-5.5*a, 5.5*a; -3.5*h, 3.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
	slab_yn = Box([-5.5*a, 5.5*a; -3.5*h, -0.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
	slab_yp = Box([-5.5*a, 5.5*a; 0.5*h, 3.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);

	hole = CircularCylinder(Axis.z, [0 0 0], r, t, [2*r/dd, 2*r/dd, t/td]);

	hole_yn_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yn);
	hole_yp_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yp);
	hole_array = [hole_yn_array, hole_yp_array];

	%% Solve the system.
	gray = [0.5 0.5 0.5];  % [r g  b]
	solveropts.method = 'gpu';

% 	[E, H, obj_array, err] = maxwell_run(c, 1e-9, 1550, ...
	[E, H, obj_array, err] = maxwell_run(1e-9, 1550, ...
		{'vacuum', 'white', 1.0}, [-5.5*a, 5.5*a; -3.5*h, 3.5*h; -3*t, 3*t], [a/ad, a/ad, t/td], BC.p, [2*a 0 t], ...
		{'Palik/Si', gray}, slab, ...
		{'vacuum', 'white', 1.0}, periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yn), ...
		{'vacuum', 'white', 1.0}, periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yp), ...
		PointSrc(Axis.y, [0, 0, 0]), solveropts, inspect_only);

	save('EH_PhC_basic.mat', 'E', 'H', 'obj_array');
else
	load('EH_PhC_basic.mat');
end

%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts
	opts.cscale = 5e-3;
% 	opts.withobj = true;
	opts.withobj = false;
%  	opts.withgrid = true;
% 	opts.withinterp = false;
	visall(E{Axis.y}, obj_array, opts);

%%
	Sx = poynting(Axis.x, E{Axis.y}, E{Axis.z}, H{Axis.y}, H{Axis.z}, Axis.x, 100);
	figure;
	vis2d(Sx, obj_array);

	p = flux_patch(Sx)
end
