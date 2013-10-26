clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Solve the system.
s = 220;
w = 20;
dL = 5;
dl = 1;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, 200, ...
	'DOM', {'vacuum', 'none', 1.0}, [-310 310; -310 310; 0 dl], [dL dL dl], BC.p, [10*dL 10*dL 0], ...
	'OBJ', ...
		{'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 dl], dl), ...
		{'Johnson/Au', 'y'}, ...
			PolygonalCylinder(Axis.z, dl, dl/2, [w/2 0; w/2+s*sqrt(3)/2 -s/2; w/2+s*sqrt(3)/2 s/2], dl), ...
			PolygonalCylinder(Axis.z, dl, dl/2, [-w/2 0; -w/2-s*sqrt(3)/2 s/2; -w/2-s*sqrt(3)/2 -s/2], dl), ...
	'SRCJ', PointSrc(Axis.x, [0 0 0]), ...
	inspect_only);

%% Visualize the solution.
figure
clear opts
opts.withinterp = false;
opts.withobjsrc = true;
opts.withabs = true;
% opts.withgrid = true;
opts.cscale = 1e-3;
z_location = 0;
vis2d(E{Axis.x}, Axis.z, z_location, obj_array, src_array, opts)
% vis2d(H{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)

% %% Calculate the power emanating from the source.
% power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
% fprintf('power = %e\n', power);
