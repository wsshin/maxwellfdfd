clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;  % true to inspect structure before any calculation

%% Solve the system.
L0 = 1e-9;  % unit of length: nm
dL = 10;  % coarse grid cell size
dl = 2;  % fine grid cell size

r = 100;  % radius
s = 200;  % semiside of TF/SF box
b = s + 200 + 10*dL;  % semiside of simulation domain (10 grid cells inside PML)
sm = s + 100;

wvlen = 1550;  % wavelength

[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', L0, wvlen, ...
	'DOM', {'vacuum', 'none', 1.0}, [-b b; -b b; 0 dL], dL, BC.p, [10*dL 10*dL 0], ...
	'OBJ', {'vacuum', 'blue', 1.0}, Box([-sm sm; -sm sm; 0 dL]), ...
	'SOBJ', ...  % scattering objects
		{'Johnson/Au', 'y'}, ...
			CircularCylinder(Axis.z, dL, [0 0 dL/2], r, [dl dl dL]), ...
	'SRCJ', TFSFPlaneSrc([-s s; -s s; 0 dL], Axis.y, Axis.x), ...  % normal direction is y; polarization direction is x
	inspect_only);


if ~inspect_only
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

	%% Calculate the scatterd power
	s_power = powerflux_box(E,H, [-sm sm; -sm sm; 0 dL]);
	fprintf('scattering cross section = %e\n', s_power/dL);  % intensity of incident wave = 0.5
end
