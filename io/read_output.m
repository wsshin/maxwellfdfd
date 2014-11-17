function [E_cell, H_cell] = read_output(filenamebase)

chkarg(istypesizeof(filenamebase, 'char', [1 0]), '"filenamebase" should be string.');

%% Read the grid.
[grid3d, ge, osc] = read_grid(filenamebase);

%% Read the solution files.
efile = ['./', filenamebase, '.E.h5'];
hfile = ['./', filenamebase, '.H.h5'];

E_array = h5read(efile, '/E');  % E_array(ri, w, i, j, k)
H_array = h5read(hfile, '/H');  % H_array(ri, w, i, j, k)
assert(all(size(E_array) == size(H_array)));

E_array = collapse_complex(E_array);  % E_array(w, i, j, k) with leading singleton dimension
H_array = collapse_complex(H_array);  % H_array(w, i, j, k) with leading singleton dimension

% Below, shiftdim() is to remove the leading singleton dimension.
E = {shiftdim(E_array(Axis.x,:,:,:)), shiftdim(E_array(Axis.y,:,:,:)), shiftdim(E_array(Axis.z,:,:,:))};
H = {shiftdim(H_array(Axis.x,:,:,:)), shiftdim(H_array(Axis.y,:,:,:)), shiftdim(H_array(Axis.z,:,:,:))};

E_cell = cell(1, Axis.count);
H_cell = cell(1, Axis.count);
for w = Axis.elems
	E_cell{w} = array2scalar(E{w}, grid3d, ge, FT.e, w, osc, PhysQ.E);
	H_cell{w} = array2scalar(H{w}, grid3d, ge, FT.h, w, osc, PhysQ.H);
end
