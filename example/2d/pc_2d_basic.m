clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%%
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

rod = CircularCylinder(Axis.z, t, [0 0 t/2], r, [2*r/dd, 2*r/dd, t]);

%% Solve the system.
gray = [0.5 0.5 0.5];  % [r g b]
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, 1550, ...
	'DOM', {'vacuum', 'white', 1.0}, [-mx*a mx*a; -my*a my*a; 0 t], [a/ad a/ad t], BC.p, [5*a 2*a 0], ...
	'OBJ', ...
		{'dielectric', gray, 11.56}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yn), ...
		{'dielectric', gray, 11.56}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yp), ...
	'SRC', PointSrc(Axis.z, [0, 0, 0]), ...
	inspect_only);

%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts
% 	opts.withgrid = true;
% 	opts.withobjsrc = true;
	opts.withabs = false;
	opts.withpml = false;
	opts.phase = pi/2;
	vis2d(E{Axis.z}, Axis.z, 0, obj_array, src_array, opts)
end
