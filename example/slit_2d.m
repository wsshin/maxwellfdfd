clear all; close all; clear classes; clc;

%% Solve the system.
inspect_only = false;
isnew = false;

if isnew
	gray = [0.5 0.5 0.5];  % [r g b]

	% % (old code)
	% [E, H, obj_array] = maxwell_run(1e-8, 118, ...
	% 	{'vacuum', 'none', 1.0}, [-107, 107; -250, 250; 0, 1], 1, [BC.p BC.Et0 BC.p], [10 10 0],...
	% 	{['Hagemann', filesep, 'Ag'], gray}, [Box([-107, -8; -150, -50; 0, 1], 1), Box([8, 107; -150, -50; 0, 1], 1)], ...  % OBJ2
	% 	PlaneSrc(Axis.y, -200, Axis.x), inspect_only);

	% (new code)
	% [E, H, obj_array] = maxwell_run(1e-9, 1180, ...
	% 	{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, [BC.p BC.Et0 BC.p], [100 100 0],...
	% 	{['Hagemann', filesep, 'Ag'], gray}, [Box([-1070, -80; -1500, -500; 0, 10], 10), Box([80, 1070; -1500, -500; 0, 10], 10)], ...  % OBJ2
	% 	PlaneSrc(Axis.y, -2000, Axis.x), inspect_only);

	% (true PEC)
	[E, H, obj_array] = maxwell_run(1e-9, 1180, ...
		{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, [BC.p BC.Et0 BC.p], [100 100 0],...
		{'PEC', gray, inf}, [Box([-1070, -80; -1500, -500; 0, 10], 10), Box([80, 1070; -1500, -500; 0, 10], 10)], ...  % OBJ2
		PlaneSrc(Axis.y, -2000, Axis.x), inspect_only);

	% % (no metal)
	% [E, H, obj_array] = maxwell_run(1e-9, 1180, ...
	% 	{'vacuum', 'none', 1.0}, [-1070, 1070; -2500, 2500; 0, 10], 10, [BC.p BC.Et0 BC.p], [100 100 0],...
	% 	PlaneSrc(Axis.y, -2000, Axis.x), inspect_only);
	
	save('slit_2d_out', 'E', 'H', 'obj_array');
else
	load('slit_2d_out');
end

%% Visualize the solution.
figure
clear opts
opts.withabs = false;  % true: abs(solution), false: real(solution)
opts.withpml = false;  % true: show PML, false: do not show PML
z_location = 0;
vis2d(E{Axis.x}, Axis.z, z_location, obj_array, opts)

%% Calculate the power flux through the slit.
y_location = 0;
power = powerflux_patch(E, H, Axis.y, y_location);
fprintf('power = %e\n', power);

