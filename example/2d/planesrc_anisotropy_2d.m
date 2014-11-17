clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

flux_loc = 20;
wvlen = 20;
lsrc = 50/2;

%% Solve the system with Jz.
% (oblique incidence)
polarization = Axis.z;
prop = Axis.y;
angle = pi/3;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-50 50; -60 60; 0 1], 1, BC.p, [0 10 0], ...
	'OBJ', {'anisotropic material', 'c', [3 3 1]}, Box([-50 50; -60 0; 0 1]), ...
	'SRCJ', PlaneSrc(prop, lsrc, polarization, 1, angle, wvlen), ...  % PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)
	inspect_only);  

figure
clear opts;
opts.withobjsrc = true;
opts.withinterp = true;
opts.withpml = true;
vis2d(E{polarization}, Axis.z, 0, obj_array, src_array, opts)

%% Solve the system with Jz, but now for anisotropic mu.
% (oblique incidence)
polarization = Axis.z;
prop = Axis.y;
angle = pi/3;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-50 50; -60 60; 0 1], 1, BC.p, [0 10 0], ...
	'OBJ', {'anisotropic material', 'c', 1, [3 3 1]}, Box([-50 50; -60 0; 0 1]), ...
	'SRCJ', PlaneSrc(prop, lsrc, polarization, 1, angle, wvlen), ...  % PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)
	inspect_only);  

figure
clear opts;
opts.withobjsrc = true;
opts.withinterp = true;
opts.withpml = true;
vis2d(E{polarization}, Axis.z, 0, obj_array, src_array, opts)

%% Solve the system with Mz.
% (oblique incidence)
polarization = Axis.z;
prop = Axis.y;
angle = pi/3;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-50 50; -60 60; 0 1], 1, BC.p, [0 10 0], ...
	'OBJ', {'anisotropic material', 'c', [3 3 1]}, Box([-50 50; -60 0; 0 1]), ...
	'SRCM', PlaneSrc(prop, lsrc, polarization, 1, angle, wvlen), ...  % PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)
	inspect_only);  

figure
clear opts;
opts.withobjsrc = true;
opts.withinterp = true;
opts.withpml = true;
vis2d(H{polarization}, Axis.z, 0, obj_array, src_array, opts)
