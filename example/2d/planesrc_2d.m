clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

flux_loc = 20;
wvlen = 20;

%% Solve the system.
% % (incidence in x-direction)
% polarization = Axis.z;
% prop = Axis.x;
% angle = 0;
% [E, H, obj_array, src_array, J] = maxwell_run(...
% 	'OSC', 1e-9, wvlen, ...
% 	'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -50, 50; 0, 1], 1, BC.p, [10 0 0], ...
% 	'OBJ', {'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
% 	'SRCJ', PlaneSrc(prop, 0, polarization), ...  % PlaneSrc(plane_normal_axis, intercept, polarization_axis)
% 	inspect_only);

% % (incidence in y-direction)
% polarization = Axis.z;
% prop = Axis.y;
% angle = 0;
% [E, H, obj_array, src_array, J] = maxwell_run(...
% 	'OSC', 1e-9, wvlen, ...
% 	'DOM', {'vacuum', 'none', 1.0}, [-50, 50; -60, 60; 0, 1], 1, BC.p, [0 10 0], ...
% 	'OBJ', {'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
% 	'SRCJ', PlaneSrc(prop, 0, polarization), ...  % PlaneSrc(plane_normal_axis, intercept, polarization_axis)
% 	inspect_only);  

% (oblique incidence)
polarization = Axis.z;
prop = Axis.y;
angle = pi/3;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-50, 50; -60, 60; 0, 1], 1, BC.p, [0 10 0], ...
	'OBJ', {'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
	'SRCJ', PlaneSrc(prop, 0, polarization, 1, angle, wvlen), ...  % PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)
	inspect_only);  

%% Visualize the solution.
figure
clear opts;
opts.withobjsrc = true;
opts.withinterp = true;
opts.withpml = true;
vis2d(E{polarization}, Axis.z, 0, obj_array, src_array, opts)

%% Calculate the power emanating from the source and test the performance of PML.
power_measured = powerflux_patch(E, H, prop, flux_loc);
power_expected = (1/2) * (0.5 * 0.5 * (100 * 1)) / cos(angle);  % (1/2) E x H x area / cos(angle)
fprintf('power:\n');
fprintf('measured = %s\n', num2str(power_measured));
fprintf('expected = %s\n', num2str(power_expected));
fprintf('error = %s%%\n',num2str((power_measured-power_expected)/power_expected*100));
