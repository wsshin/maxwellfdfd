clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Solve the system.
flux_loc = 10;
wvlen = 20;
% % (incidence in x-direction)
% polarization = Axis.z;
% [E, H, obj_array, src_array, J] = maxwell_run(...
% 	'OSC', 1e-9, wvlen, ...
% 	'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -50, 50; 0, 1], 1, BC.p, [10 0 0], ...
% 	'OBJ', {'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
% 	'SRC', PlaneSrc(Axis.x, 0, polarization), ...  % PlaneSrc(plane_normal_axis, intercept, polarization_axis)
% 	inspect_only);

% % (incidence in y-direction)
% polarization = Axis.z;
% [E, H, obj_array, src_array, J] = maxwell_run(...
% 	'OSC', 1e-9, wvlen, ...
% 	'DOM', {'vacuum', 'none', 1.0}, [-50, 50; -60, 60; 0, 1], 1, BC.p, [0 10 0], ...
% 	'OBJ', {'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
% 	'SRC', PlaneSrc(Axis.y, 0, polarization), ...  % PlaneSrc(plane_normal_axis, intercept, polarization_axis)
% 	inspect_only);  

% (oblique incidence)
polarization = Axis.z;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-50, 50; -60, 60; 0, 1], 1, BC.p, [0 10 0], ...
	'OBJ', {'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
	'SRC', PlaneSrc(Axis.y, 0, polarization, 1, pi/4, wvlen), ...  % PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)
	inspect_only);  

%% Visualize the solution.
figure
clear opts;
opts.withobjsrc = true;
opts.withinterp = true;
vis2d(E{polarization}, Axis.z, 0, obj_array, src_array, opts)

%% Calculate the power emanating from the source.
power = powerflux_patch(E, H, Axis.y, flux_loc);
fprintf('power = %e\n', power);
