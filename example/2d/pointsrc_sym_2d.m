clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Solve the system.
% PEC symmetry plane at x = 0
[E, H, obj_array, src_array] = maxwell_run(...
	'OSC', 1e-9, 100, ...
	'DOM', {'vacuum', 'none', 1.0}, [0, 60; -60, 60; 0, 1], 1, [BC.e BC.m BC.p], [0 10; 10 10; 0 0], ...
	'SRCJ', PointSrc(Axis.z, [10, 0, 0.5]), ...
	inspect_only);
F = E;

% % PMC symmetry plane at x = 0
% [E, H, obj_array, src_array] = maxwell_run(...
% 	'OSC', 1e-9, 100, ...
% 	'DOM', {'vacuum', 'none', 1.0}, [0, 60; -60, 60; 0, 1], 1, [BC.m BC.m BC.p], [0 10; 10 10; 0 0], ...
% 	'SRCM', PointSrc(Axis.z, [10, 0, 0]), ...
% 	inspect_only);
% F = H;

%% Visualize the solution.
figure
clear opts
opts.withpml = false;
opts.withabs = true;
opts.withobjsrc = false;
vis2d(F{Axis.z}, Axis.z, 0, obj_array, src_array, opts)

% %% Calculate the power emanating from the source.
% power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
% fprintf('power = %e\n', power);
