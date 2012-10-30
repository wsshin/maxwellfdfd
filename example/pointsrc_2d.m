clear all; close all; clear classes; clc;

%% Set flags.
isnew = true;
inspect_only = false;

%% Solve the system.
if isnew
	[E, H, obj_array] = maxwell_run(1e-9, 100, ...
		{'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, BC.p, [10 10 0], ...
		PointSrc(Axis.z, [0, 0, 0]), inspect_only);

	if ~inspect_only
		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end

%% Visualize the solution.
figure
vis2d(E{Axis.z}, Axis.z, 0, obj_array)

%% Calculate the power emanating from the source.
power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
fprintf('power = %e\n', power);
