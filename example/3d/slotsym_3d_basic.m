clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = true;

%% Solve the system.
gray = [0.5 0.5 0.5];  % [r g b]
solveropts.method = 'aws';
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, 1550, ...
	'DOM', {'Palik/SiO2', 'none'}, [0, 700; 0, 600; -200, 1700], 20, [BC.e, BC.m, BC.p], [0 200; 0 200; 200 200], ...
	'OBJ', ...
		{'Palik/SiO2', 'none'}, Box([0, 50; 0, 50; -200, 1700], [2, 2, 20]), ...
		{'CRC/Ag', gray}, Box([25, 700; 0, 25; -200, 1700], 20), ...
	'SRC', ModalSrc(Axis.z, 200, 2.0), ...
	solveropts, inspect_only);

%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts;
	opts.withabs = true;
	opts.withgrid = true;
	visall(E{Axis.x}, obj_array, src_array, opts);
end
