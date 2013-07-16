clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Solve the system.
s = 120;
w = 40;
b = (s + w/2)*2;
dL = 10;
dl = 2;
wvlen = 200;
clear solveropts;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-610 610; -610 610; 0 dl], [dL dL dl], BC.p, [10*dL 10*dL 0], ...
	'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 dl], dl), ...
	'SOBJ', ...  % scatter objects
		{'Johnson/Au', 'y'}, ...
			PolygonalCylinder(Axis.z, dl, dl/2, [w/2 0; w/2+s*sqrt(3)/2 -s/2; w/2+s*sqrt(3)/2 s/2], dl), ...
			PolygonalCylinder(Axis.z, dl, dl/2, [-w/2 0; -w/2-s*sqrt(3)/2 s/2; -w/2-s*sqrt(3)/2 -s/2], dl), ...
	'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 dl], Axis.y, Axis.x), ...
	inspect_only);

%% Visualize the solution.
figure
clear opts
opts.withobjsrc = true;
opts.withabs = true;
% opts.withinterp = false;
% opts.withgrid = true;
% opts.cscale = 1e-1;
% opts.cmax = 1.4;
z_location = 0;
vis2d(E{Axis.x}, Axis.z, z_location, obj_array, src_array, opts)
% vis2d(H{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)

% %% Calculate the power emanating from the source.
% power = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
% fprintf('power = %e\n', power);
