clear all; close all; clear classes; clc;

%% Set flags.
is_new = true;
inspect_only = true;
filenamebase = 'slot_3d_basic';

if is_new
	%% Create input files.
	gray = [0.5 0.5 0.5];  % [r g b]
	solveropts.method = 'inputfile';
	solveropts.filenamebase = filenamebase;
	solveropts.showstruct = false;

	% (along x)
	modeopts.clue = 'guess';
	modeopts.neff = 2;
	[E, H, obj_array, src_array, J] = maxwell_run(...
		'OSC', 1e-9, 1550, ...
		'DOM', {'Palik/SiO2', 'none'}, [-200 1700; -700 700; -600 600], 20, BC.p, 200, ...
		'OBJ', ...
			{'Palik/SiO2', 'none'}, Box([-200 1700; -50 50; -50 50], [20 2 2]), ...
			{'CRC/Ag', gray}, ...
				Box([-200 1700; -700 -25; -25 25], 20), Box([-200 1700; 25 700; -25 25], 20), ...
		'SRCJ', ModalSrc(Axis.x, 200, modeopts), ...
		solveropts, inspect_only);
	
	save(filenamebase, 'obj_array', 'src_array');
else
	[E, H] = read_output(filenamebase);
	load(filenamebase);
end

%% Visualize the solution.
if ~is_new
	figure;
	clear opts;
	opts.withabs = true;
	opts.isopaque = false;
% 	opts.withgrid = true;
	visall(E{Axis.x}, obj_array, src_array, opts);
end
