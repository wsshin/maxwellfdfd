clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = true;

%% Create shapes.
a = 420;  % lattice constant
t = 0.6*a;  % slab thickness
r = 0.29*a;  % hole radius
h = sqrt(3)/2*a;  % distance between rows of holes

ad = 25;  % divider for a
td = 10;  % divider for t
dd = 10;  % divider for d = 2*r

slab = Box([-5.5*a, 5.5*a; -3.5*h, 3.5*h; -t/2, t/2]);
slab_yn = Box([-5.5*a, 5.5*a; -3.5*h, -0.5*h; -t/2, t/2]);
slab_yp = Box([-5.5*a, 5.5*a; 0.5*h, 3.5*h; -t/2, t/2]);

hole = CircularCylinder(Axis.z, t, [0 0 0], r);

hole_yn_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yn);
hole_yp_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yp);
hole_array = [hole_yn_array, hole_yp_array];

%% Solve the system.
gray = [0.5 0.5 0.5];  % [r g b]
solveropts.method = 'aws';
withuniformgrid = true;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, 1550, ...
	'DOM', {'vacuum', 'white', 1.0}, [-5.5*a, 5.5*a; -3.5*h, 3.5*h; -3*t, 3*t], [11*a/220, 7*h/71, t/td], BC.p, [2*a 0 t], withuniformgrid, ...
	'OBJ', ...
		{'Palik/Si', gray}, slab, ...
		{'vacuum', 'white', 1.0}, periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yn), ...
		{'vacuum', 'white', 1.0}, periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yp), ...
	'SRCJ', PointSrc(Axis.y, [0, 0, 0]), ...
	solveropts, inspect_only);

%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts;
	opts.cscale = 5e-3;
	opts.withobjsrc = false;
	visall(E{Axis.y}, obj_array, src_array, opts);
end
