clear all; close all; clear classes; clc;

%% Set flags.
isnew = true;
inspect_only = false;

%% Solve the system.
if isnew
	[E, H, obj_array, src_array] = maxwell_run(...
		'OSC', 1e-9, 100, ...
		'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, BC.p, [10 10 0], ...
		'SRC', PointSrc(Axis.z, [0, 0, 0]), ...
		inspect_only);
	z_location = 0;

% 	[E, H, obj_array, src_array] = maxwell_run(...
% 		'OSC', 1e-9, 100, ...
% 		'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, BC.p, [10 10 0], ...
% 		'SRC', PointSrcM(Axis.z, [0, 0, 0.5]), ...
% 		inspect_only);
% 	z_location = 0.5;

	if ~inspect_only
		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end

%% Visualize the solution.
figure
clear opts
opts.withinterp = false;
opts.withobjsrc = true;
% opts.cscale = 1e-3;
vis2d(E{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)
% vis2d(H{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)

% %% Calculate the power emanating from the source.
% power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
% fprintf('power = %e\n', power);
