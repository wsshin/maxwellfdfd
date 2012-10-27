clear all; close all; clear classes; clc;
cluster = 'wsshin2';
nodes = 2;

%% Initialize.
maxwell.aws_credentials('AKIAJJUWULDKPRQFDXFA', 'otTKcW/ivRjTkvaV9ZlZBNALpzkABU1cSI6vTKDq')
maxwell.launch(cluster, nodes);

%% Solve.
gray = [0.5 0.5 0.5];  % [r g b]
inspect_only = false;
solveropts.method = 'aws';
solveropts.cluster = cluster;
solveropts.nodes = nodes;
[E, H, obj_array] = maxwell_run(1e-9, 1550, ...
		{['Palik', filesep, 'SiO2'], 'none'}, [-700, 700; -600, 600; -200, 1700], 20, BC.p, 200, ...
		{['Palik', filesep, 'SiO2'], 'none'}, Box([-50, 50; -50, 50; -200, 1700], [2, 2, 20]), ...
		{['Hagemann', filesep, 'Ag'], gray}, ...
		[Box([-700, -25; -25, 25; -200, 1700], 20), Box([25, 700; -25, 25; -200, 1700], 20)], ...
		DistributedSrc(Axis.z, 200, 2.0), 1e-5, solveropts, inspect_only);

%% Terminate.
maxwell.terminate(cluster)

%% Visualize.
if ~inspect_only
	figure;
	clear opts
	opts.withabs = true;
	opts.withobj = true;
	visall(E{Axis.x}, obj_array, opts);
end
