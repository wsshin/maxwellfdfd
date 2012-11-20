clear all; close all; clear classes; clc;

%% Set flags.
isnew = true;
inspect_only = false;

%% Solve the system.
flux_loc = 10;
if isnew
	[E, H, obj_array, src_array] = maxwell_run(1e-9, 100, ...
		{'vacuum', 'none', 1.0}, [-50, 50; -60, 60; 0, 1], 1, BC.p, [0 10 0], ...
		{'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
		PlaneSrc(Axis.y, 0, Axis.z), inspect_only);  % PlaneSrc(plane_normal_axis, intercept, polarization_axis)

% 	[E, H, obj_array, src_array] = maxwell_run(1e-9, 100, ...
% 		{'vacuum', 'none', 1.0}, [-60, 60; -50, 50; 0, 1], 1, BC.p, [10 0 0], ...
% 		{'vacuum', 'none', 1.0}, Plane(Axis.y, flux_loc), ...
% 		PlaneSrc(Axis.x, 0, Axis.z), inspect_only);  % PlaneSrc(plane_normal_axis, intercept, polarization_axis)

% 	[E, H, obj_array, src_array] = maxwell_run(1e-9, 100, ...
% 		{'vacuum', 'none', 1.0}, [0, 1; -50, 50; -60, 60], 1, BC.p, [0 0 10], ...
% 		{'vacuum', 'none', 1.0}, Plane(Axis.z, flux_loc), ...
% 		PlaneSrc(Axis.z, 0, Axis.x), inspect_only);  % PlaneSrc(plane_normal_axis, intercept, polarization_axis)

	
	if ~inspect_only
		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end

%% Visualize the solution.
figure
clear opts;
opts.withobj = false;
vis2d(E{Axis.z}, Axis.z, 0, obj_array, src_array, opts)
% vis2d(E{Axis.x}, Axis.x, 0, obj_array, src_array)

%% Calculate the power emanating from the source.
power = powerflux_patch(E, H, Axis.y, flux_loc);
fprintf('power = %e\n', power);
