clear all; close all; clear classes; clc;

%% Construct parameters.
L0 = 1e-9;
wvlen = 1550;
unit = PhysUnit(L0);
osc = Oscillation(wvlen, unit);

silica = Material.create(['Palik', filesep, 'SiO2'], 'none', osc);
Ag = Material.create(['Hagemann', filesep, 'Ag'], [0.5 0.5 0.5], osc);  % [0.5 0.5 0.5]: gray in RGB

domain = Domain([-700, 700; -600, 600; -200, 1700], 20);
domain_silica = Object(domain, silica);

refined_domain = Box([-50, 50; -50, 50; -200, 1700], [2, 2, 20]);
refined_domain_silica = Object(refined_domain, silica);

film1 = Box([-700, -25; -25, 25; -200, 1700], 20);
film1_Ag = Object(film1, Ag);

film2 = Box([25, 700; -25, 25; -200, 1700], 20);
film2_Ag = Object(film2, Ag);

src = DistributedSrc(Axis.z, 200, 2.0);


%% Solve the system.
inspect_only = true;
[E, H, obj_array, err] = ...
	maxwell_run(osc, domain_silica, [BC.p BC.p; BC.p BC.p; BC.p BC.p], [200 200 200], refined_domain_silica, film1_Ag, film2_Ag, src, inspect_only);


%% Visualize the solution.
if ~inspect_only
	figure;
	opts.withabs = true;
	visall(E{Axis.x}, obj_array, opts);
end
