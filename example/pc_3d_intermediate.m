clear all; close all; clear classes; clc;

%% Set flags.
isnew = true;
inspect_only = false;

%%
L0 = 1e-9;
wvlen = 1550;
unit = PhysUnit(L0);
osc = Oscillation(wvlen, unit);

%%
vac = Material('vacuum', 'w', 1.0);
Si = Material.create(['Palik', filesep, 'Si'], [0.5 0.5 0.5], osc);  % [0.5 0.5 0.5]: gray in RGB

%% Create shapes.
a = 420;  % lattice constant
t = 0.6*a;  % slab thickness
r = 0.29*a;  % hole radius
h = sqrt(3)/2*a;  % distance between rows of holes

ad = 25;  % divider for a
td = 10;  % divider for t
dd = 10;  % divider for d = 2*r

domain = Domain([-5.5*a, 5.5*a; -3.5*h, 3.5*h; -3*t, 3*t], [a/ad, a/ad, t/td]);
domain_vac = Object(domain, vac);

slab = Box([-5.5*a, 5.5*a; -3.5*h, 3.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
slab_Si = Object(slab, Si);

hole = CircularCylinder(Axis.z, [0 0 0], r, t, [2*r/dd, 2*r/dd, t/td]);

slab_yn = Box([-5.5*a, 5.5*a; -3.5*h, -0.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
hole_yn_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yn);

slab_yp = Box([-5.5*a, 5.5*a; 0.5*h, 3.5*h; -t/2, t/2], [a/ad, a/ad, t/td]);
hole_yp_array = periodize_shape(hole, {[a 0 0], [a/2 h 0], [0 0 t]}, slab_yp);

hole_array = [hole_yn_array, hole_yp_array];
hole_array_vac = Object(hole_array, vac);

src = PointSrc(Axis.y, [0, 0, 0]);

%% Solve the system.
if isnew
	gray = [0.5 0.5 0.5];  % [r g b]
	[E, H, obj_array, err] = maxwell_run(...
		'OSC', osc, ...
		'DOM', domain_vac, BC.p, [2*a 0 t], ...
		'OBJ', slab_Si, hole_array_vac, ...
		'SRC', src, ...
		inspect_only);

	if ~inspect_only
		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end

%% Visualize the solution.
if ~inspect_only
	figure;
	opts.cscale = 5e-3;
	visall(E{Axis.y}, obj_array, opts);
end
