clear all; close all; clear classes; clc;

%% Set flags.
isnew = true;
inspect_only = false;

%% Solve the system.
if isnew
	gray = [0.5 0.5 0.5];  % [r g b]
	flux_y = -1000;
	flux_x1 = -200; flux_x2 = 200;

% 	% (old code)
% 	[E, H, obj_array, src_array] = maxwell_run(1e-8, 118, ...
% 		{'vacuum', 'none', 1.0}, [-107, 107; -250, 250; 0, 1], 1, BC.p, [10 10 0],...
% 		{'vacuum', 'b', 1.0}, Rectangle(Axis.y, flux_y, [0 10; flux_x1 flux_x2]), ...
% 		{['CRC', filesep, 'Ag'], gray}, Box([-107, -8; -150, -50; 0, 1]), Box([8, 107; -150, -50; 0, 1]), ...  % metal slit
% 		PlaneSrc(Axis.y, -200, Axis.x), inspect_only);

% 	% (new code)
% 	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1180, ...
% 		{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, BC.p, [100 100 0],...
% 		{'vacuum', 'b', 1.0}, Rectangle(Axis.y, flux_y, [0 10; flux_x1 flux_x2]), ...
% 		{['CRC', filesep, 'Ag'], gray}, Box([-1070, -80; -1500, -500; 0, 10]), Box([80, 1070; -1500, -500; 0, 10]), ...  % metal slit
% 		PlaneSrc(Axis.y, -2000, Axis.x), inspect_only);

% 	% (true PEC)
% 	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1180, ...
% 		{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, BC.p, [100 100 0],...
% 		{'vacuum', 'b', 1.0}, Rectangle(Axis.y, flux_y, [0 10; flux_x1 flux_x2]), ...
% 		{'PEC', gray, Inf}, Box([-1070, -80; -1500, -500; 0, 10]), Box([80, 1070; -1500, -500; 0, 10]), ...  % metal slit
% 		PlaneSrc(Axis.y, -2000, Axis.x), inspect_only);

% 	% (no metal)
% 	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1180, ...
% 		{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, BC.p, [100 100 0],...
% 		{'vacuum', 'b', 1.0}, Plane(Axis.y, flux_y), ...
% 		PlaneSrc(Axis.y, -2000, Axis.x), inspect_only);

% 	% (distributed source)
% 	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1180, ...
% 		{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, BC.p, [100 100 0],...
% 		{'vacuum', 'b', 1.0}, Rectangle(Axis.y, flux_y, [0 10; flux_x1 flux_x2]), ...
% 		{['CRC', filesep, 'Ag'], gray}, Box([-1070, -80; -1500, -500; 0, 10]), Box([80, 1070; -1500, -500; 0, 10]), ...  % metal slit
% 		DistributedSrc(Axis.y, -1000, 1.0), inspect_only);

% 	% (point source)
% 	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1180, ...
% 		{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, BC.p, [100 100 0],...
% 		{'vacuum', 'b', 1.0}, Rectangle(Axis.y, flux_y, [0 10; flux_x1 flux_x2]), ...
% 		{['CRC', filesep, 'Ag'], gray}, Box([-1070, -80; -1500, -500; 0, 10]), Box([80, 1070; -1500, -500; 0, 10]), ...  % metal slit
% 		PointSrc(Axis.x, [0, -2000, 5]), inspect_only);

	% (rectangular source)
	xc = 101;
	[E, H, obj_array, src_array] = maxwell_run(1e-9, 1180, ...
		{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, BC.p, [100 100 0],...
		{'vacuum', 'b', 1.0}, Rectangle(Axis.y, flux_y, [0 10; flux_x1 flux_x2]), ...
		{['CRC', filesep, 'Ag'], gray}, Box([-1070, -80; -1500, -500; 0, 10]), Box([80, 1070; -1500, -500; 0, 10]), ...  % metal slit
		RectSrc(Axis.y, -2000, [0 10; xc-80 xc+80], Axis.x), inspect_only);
	if ~inspect_only
		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end

if ~inspect_only
	%% Visualize the solution.
	figure
	clear opts
	opts.withabs = true;  % true: abs(solution), false: real(solution)
	opts.withpml = false;  % true: show PML, false: do not show PML
	opts.withgrid = false;
	opts.cscale = 1e-1;
	z_loc = 0;
	vis2d(E{Axis.x}, Axis.z, z_loc, obj_array, src_array, opts)

	%% Calculate the power flux through the slit.
	power = powerflux_patch(E, H, Axis.y, flux_y, [0 10; flux_x1 flux_x2]);
	fprintf('power = %e\n', power);
end
