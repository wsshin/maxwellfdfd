clear all; close all; clear classes; 

%% Set flags.
isnew = true;
inspect_only = false;

%% Solve the system.
if isnew
	s = 120;
	w = 40;
	b = (s + w/2)*2;
	dL = 10;
	dl = 2;
	wvlen = 200;
	solveropts.withinterp = false;
	[E, H, obj_array, src_array, J] = maxwell_run(...
		'OSC', 1e-9, wvlen, ...
		'DOM', {'vacuum', 'none', 1.0}, [-610 610; -610 610; 0 1], [dL dL 1], BC.p, [10*dL 10*dL 0], ...
		'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], [dl dl 1]), ...
		'SOBJ', ...  % scatter objects
			{['CRC', filesep, 'Ag'], 'y'}, ...
				PolygonalCylinder(Axis.z, 1, 1/2, [w/2 0; w/2+s*sqrt(3)/2 -s/2; w/2+s*sqrt(3)/2 s/2], [dl dl 1]), ...
				PolygonalCylinder(Axis.z, 1, 1/2, [-w/2 0; -w/2-s*sqrt(3)/2 s/2; -w/2-s*sqrt(3)/2 -s/2], [dl dl 1]), ...
		'SRC', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
		solveropts, inspect_only);

	if ~inspect_only
% 		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end

%% Visualize the solution.
figure
clear opts
opts.withobjsrc = true;
opts.withabs = true;
% opts.withinterp = false;
% opts.withgrid = true;
% opts.cscale = 1e-1;
z_location = 0;
% vis2d(E{Axis.x}, Axis.z, z_location, obj_array, src_array, opts)
% vis2d(H{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)

% %% Calculate the power emanating from the source.
% Exmax = max(abs(E{Axis.x}.array(:)));
% Hzmax = max(abs(H{Axis.z}.array(:)));
Exmax = max(abs(E{Axis.x}(:)));
Hzmax = max(abs(H{Axis.z}(:)));
fprintf('Exmax = %e, Hzmax = %e\n', Exmax, Hzmax);
% power = powerflux_box(E,H,[-b-10 b+10; -b-10 b+10; 0 1]);
% fprintf('power = %e\n', power);
