clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Solve the system.
wvlen = 20;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, BC.p, [10 10 0], ...
	'SRC', PointSrc(Axis.z, [0, 0, 0]), ...
	inspect_only);

% [E, H, obj_array, src_array] = maxwell_run(...
% 	'OSC', 1e-9, wvlen, ...
% 	'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, BC.p, [10 10 0], ...
% 	'SRC', PointSrcM(Axis.z, [0, 0, 0.5]), ...
% 	inspect_only);

%% Visualize the solution.
figure
clear opts
opts.withinterp = true;
opts.withobjsrc = true;
opts.withpml = true;
% opts.cscale = 1e-2;
vis2d(E{Axis.z}, Axis.z, 0, obj_array, src_array, opts)
% vis2d(H{Axis.z}, Axis.z, 0, obj_array, src_array, opts)

% %% Calculate the power emanating from the source.
% power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
% fprintf('power = %e\n', power);
