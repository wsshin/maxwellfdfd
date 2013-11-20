clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = true;

%% Solve the system.
gray = [0.5 0.5 0.5];  % [r g b]
solveropts.method = 'inputfile';
solveropts.filenamebase = 'slot_3d';

% (along x)
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, 1550, ...
	'DOM', {'Palik/SiO2', 'none'}, [-200 1700; -700 700; -600 600], 20, BC.p, 200, ...
	'OBJ', ...
		{'Palik/SiO2', 'none'}, Box([-200 1700; -50 50; -50 50], [20 2 2]), ...
		{'CRC/Ag', gray}, ...
			Box([-200 1700; -700 -25; -25 25], 20), Box([-200 1700; 25 700; -25 25], 20), ...
	'SRCJ', ModalSrc(Axis.x, 200), ...
	solveropts, inspect_only);

%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts;
	opts.withabs = true;
% 	opts.withgrid = true;
	visall(E{Axis.x}, obj_array, src_array, opts);
end
