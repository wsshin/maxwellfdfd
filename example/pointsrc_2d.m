clear all; close all; clear classes; clc;

%% Solve the system.
isnew = false;
inspect_only = false;

if isnew
	[E, H, obj_array] = maxwell_run(1e-9, 100, ...
		{'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, [BC.Et0 BC.Et0 BC.p], [10 10 0], ...
		PointSrc(Axis.z, [0, 0, 0.5]), inspect_only);
	save('pointsrc_2d_out', 'E', 'H', 'obj_array');
else
	load('pointsrc_2d_out');
end

%% Visualize the solution.
figure
vis2d(E{Axis.z}, Axis.z, 0, obj_array)

%% Calculate the power emanating from the source.
power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
fprintf('power = %e\n', power);
