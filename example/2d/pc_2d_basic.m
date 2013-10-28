clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Create shapes.
a = 420;  % lattice constant
t = 1;  % slab thickness

% permittivity, r/a, and omega*a/(2*pi*c) are taken from p.234 of John D.
% Joannopoulos? et al., "Photonic Crystals: Molding the Flow of Light," 2nd
% edition so that the frequency lies in a band gap.
eps_diel = 11.4;
r = 0.25*a;  % hole radius
wvlen = a/0.3;  % omega*a/(2*pi*c) = k*a/(2*pi) = a/lambda = 0.3

ad = 25;  % divider for a
td = 10;  % divider for t
dd = 10;  % divider for d = 2*r

mx = 25.5;  % half integer puts domain boundary between cylinders and makes PML works better
my = 7.5;
slab_yn = Box([-mx*a mx*a; -my*a -0.5*a; 0 t], [a/ad, a/ad, t]);
slab_yp = Box([-mx*a mx*a; 0.5*a my*a; 0 t], [a/ad, a/ad, t]);

rod = CircularCylinder(Axis.z, t, [0 0 t/2], r, [2*r/dd, 2*r/dd, t]);

%% Solve the system.
gray = [0.5 0.5 0.5];  % [r g b]
src_loc = 6*a;
[E, H, obj_array, src_array, J] = maxwell_run(...
	'OSC', 1e-9, wvlen, ...
	'DOM', {'vacuum', 'white', 1.0}, [-mx*a mx*a; -my*a my*a; 0 t], [a/ad a/ad t], BC.p, [7*a 2*a 0], 2, 1e-4, ...
	'OBJ', ...
		{'dielectric', gray, eps_diel}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yn), ...
		{'dielectric', gray, eps_diel}, periodize_shape(rod, {[a 0 0], [0 a 0], [0 0 t]}, slab_yp), ...
	'SRCJ', PointSrc(Axis.z, [src_loc, 0, 0.5]), ...
	inspect_only);

%% Visualize the solution.
if ~inspect_only
	figure;
	clear opts
% 	opts.withgrid = true;
	opts.withobjsrc = false;
	opts.withabs = true;
	opts.withpml = false;
	opts.phase = pi/2;
	figure(1)
	vis2d(E{Axis.z}, Axis.z, 0.5, obj_array, src_array, opts)
	%%
	figure(2)
	vis2d(H{Axis.y}, Axis.z, 0.5, obj_array, src_array, opts)
	
	%%
	flux_loc = 3*a;
	power_right = powerflux_patch(E, H, Axis.x, src_loc + flux_loc);
	power_left = -powerflux_patch(E, H, Axis.x, src_loc - flux_loc);
	fprintf('power:\n');
	fprintf('right = %s\n', num2str(power_right));
	fprintf('left = %s\n', num2str(power_left));
	fprintf('error = %s%%\n',num2str((power_left-power_right)/power_right*100));
	
	%%
	Sx = poynting(Axis.x, E{Axis.y}, E{Axis.z}, H{Axis.y}, H{Axis.z}, Axis.y, 0);
	[array, l] = Sx.data_original;
	figure(3)
	plot(l{2}, abs(array))
	mx*a - 10*a
	
% 	%%
% 	xs = E{Axis.z}.grid3d.l{Axis.x,GT.prim};
% 	ez = NaN(size(xs));
% 	i = 0;
% 	for x = xs
% 		i = i + 1;
% 		ez(i) = E{Axis.z}.value([x 0 t/2]);
% 	end
% 	%%
% 	figure(4);
% 	plot(xs, imag(ez));
% 
% 	%%
% 	xs = E{Axis.z}.grid3d.l{Axis.x,GT.prim};
% 	hy = NaN(size(xs));
% 	i = 0;
% 	for x = xs
% 		i = i + 1;
% 		hy(i) = H{Axis.y}.value([x 0 t/2]);
% 	end
% 	%%
% 	figure(5);
% 	plot(xs, real(hy));
% 	
% 	%%
% % 	sx = real(imag(ez) .* real(hy));
% % 	sx = imag(ez .* conj(hy));
% 	sx = real(ez .* conj(hy));
% 	figure(6);
% 	plot(xs, abs(sx));
% 	
% 	% Need to verify if the same problem occurs for precisely selected frequency
% 	% for defect waveguide.
end
