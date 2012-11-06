clear all; close all; clear classes; clc;

%% Set flags.
isnew = true;
inspect_only = false;
cluster = 'wsshin';
nodes = 2;

%%
if isnew
	%% Initialize.
	if ~inspect_only
		maxwell.aws_credentials('AKIAJJUWULDKPRQFDXFA', 'otTKcW/ivRjTkvaV9ZlZBNALpzkABU1cSI6vTKDq')
		maxwell.launch(cluster, nodes);
	end

	%% Solve.
	gray = [0.5 0.5 0.5];  % [r g b]
	solveropts.method = 'aws';
	solveropts.cluster = cluster;
	solveropts.nodes = nodes;
	[E, H, obj_array, err] = maxwell_run(1e-9, 1550, ...
		{['Palik', filesep, 'SiO2'], 'none'}, [-700, 700; -600, 600; -200, 1700], 20, BC.p, 200, ...
		{['Palik', filesep, 'SiO2'], 'none'}, Box([-50, 50; -50, 50; -200, 1700], [2, 2, 20]), ...
		{['CRC', filesep, 'Ag'], gray}, ...
		[Box([-700, -25; -25, 25; -200, 1700], 20), Box([25, 700; -25, 25; -200, 1700], 20)], ...
		DistributedSrc(Axis.z, 200, 2.0), solveropts, inspect_only);

	if ~inspect_only
		save(mfilename, 'E', 'H', 'obj_array');
	end
else
	load(mfilename);
end

%% Visualize.
if ~inspect_only
	figure;
	clear opts
	opts.withabs = false;
	opts.withobj = true;
	visall(E{Axis.x}, obj_array, opts);
end

if false
	%% Terminate.
	maxwell.terminate(cluster)
end
