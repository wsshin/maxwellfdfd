clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = true;

%%
L0 = 1e-9;
wvlen = 1550;
unit = PhysUnit(L0);
osc = Oscillation(wvlen, unit);

%% Create materials.
vac = Material('vacuum', 'w', 1.0);
Si = Material.fromtable(osc, 'Palik/Si', [0.5 0.5 0.5]);  % [0.5 0.5 0.5]: gray in RGB

%% Create objects.
a = 420;  % lattice constant
t = 0.6*a;  % slab thickness
r = 0.29*a;  % hole radius
h = sqrt(3)/2*a;  % distance between rows of holes

ad = 25;  % divider for a
td = 10;  % divider for t
dd = 10;  % divider for d = 2*r

domain = Domain([-5.5*a, 5.5*a; -3.5*h, 3.5*h; -3*t, 3*t], [a/ad, a/ad, t/td]);
domain_vac = EMObject(domain, vac);

slab = Box([-5.5*a, 5.5*a; -3.5*h, 3.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
slab_Si = EMObject(slab, Si);

hole = CircularCylinder(Axis.z, t, [0 0 0], r, [2*r/dd, 2*r/dd, t/td]);

slab_yn = Box([-5.5*a, 5.5*a; -3.5*h, -0.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
hole_yn_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yn);

slab_yp = Box([-5.5*a, 5.5*a; 0.5*h, 3.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
hole_yp_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yp);

hole_array = [hole_yn_array, hole_yp_array];
hole_array_vac = EMObject(hole_array, vac);

%% Create a source.
src = PointSrc(Axis.y, [0, 0, 0]);

%% Solve the system.
gray = [0.5 0.5 0.5];  % [r g b]
solveropts.method = 'inputfile';
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', osc, ...
	'DOM', domain_vac, BC.p, [2*a 0 t], ...
	'OBJ', slab_Si, hole_array_vac, ...
	'SRCJ', src, ...
	solveropts, inspect_only);

%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts;
	opts.cscale = 5e-3;
	visall(E{Axis.y}, obj_array, src_array, opts);
end
