clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Solve the system.
wvlen = 20;

% (electric current source)
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, BC.p, [10 10 0], ...
	'SRCJ', PointSrc(Axis.z, [0, 0, 0.5]), ...
	inspect_only);

% % (magnetic current source)
% [E, H, obj_array, src_array] = maxwell_run(...
% 	'OSC', 1e-9, wvlen, ...
% 	'DOM', {'vacuum', 'none', 1.0}, [-60, 60; -60, 60; 0, 1], 1, BC.p, [10 10 0], ...
% 	'SRCM', PointSrc(Axis.z, [0, 0, 0]), ...
% 	inspect_only);

%% Visualize the solution.
figure
clear opts
% opts.withinterp = false;
opts.withobjsrc = true;
% opts.withpml = false;
% opts.withabs = true;
opts.cscale = 2e-1;
vis2d(E{Axis.z}, Axis.z, 0.5, obj_array, src_array, opts)
% vis2d(H{Axis.x}, Axis.z, 0.5, obj_array, src_array, opts)

% %% Calculate the power emanating from the source.
% power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
% fprintf('power = %e\n', power);
