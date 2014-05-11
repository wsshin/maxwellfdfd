function [E_cell, H_cell] = read_output(filenamebase)

chkarg(istypesizeof(filenamebase, 'char', [1 0]), '"filenamebase" should be string.');

%% Read the grid.
[grid3d, ge, osc] = read_grid(filenamebase);

%% Read the solution files.
efile = ['./', filenamebase, '.E.h5'];
hfile = ['./', filenamebase, '.H.h5'];

E_array = h5read(efile, '/E');
H_array = h5read(hfile, '/H');
assert(all(size(E_array) == size(H_array)));

E_array = collapse_complex(E_array);
H_array = collapse_complex(H_array);

%% Permute E_array and H_array's indices to (x, y, z, axis)
nD = ndims(E_array);
E_array = permute(E_array, [2:nD, 1]);
H_array = permute(H_array, [2:nD, 1]);

E_cell = cell(1, Axis.count);
H_cell = cell(1, Axis.count);
ind = cell(1, nD);
for w = Axis.elems
	for n = 1:nD
		ind{n} = ':';
	end
	ind{end} = int(w);

	gt = ge;  % grid type for E-field
	gt_array = gt(ones(1, Axis.count));
	gt_array(w) = alter(gt);
	E_cell{w} = array2scalar(E_array(ind{:}), PhysQ.E, grid3d, w, FT.e, gt_array, osc);

	gt = alter(ge);  % grid type for H-field
	gt_array = gt(ones(1, Axis.count));
	gt_array(w) = alter(gt);
	H_cell{w} = array2scalar(H_array(ind{:}), PhysQ.H, grid3d, w, FT.h, gt_array, osc);
end
