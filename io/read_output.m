function [E_cell, H_cell] = read_output(filenamebase)
chkarg(istypesizeof(filenamebase, 'char', [1 0]), '"filenamebase" should be string.');

inputfile = ['./', filenamebase, '.h5'];
efile = ['./', filenamebase, '.E.h5'];
hfile = ['./', filenamebase, '.H.h5'];

lambda = h5read(inputfile, '/lambda');
L0 = h5read(inputfile, '/L0');
Npml = h5read(inputfile, '/Npml').';
bc = BC.elems(h5read(inputfile, '/bc').' + 1);
bc = bc(:, Sign.n).';

lprim = cell(1, Axis.count);
for w = Axis.elems
	lprim{w} = h5read(inputfile, ['/', char(w), '_prim']).';
end

unit = PhysUnit(L0);
osc = Oscillation(lambda, unit);
grid3d = Grid3d(unit, lprim, Npml, bc);

E_array = h5read(efile, '/E');
H_array = h5read(hfile, '/H');
assert(all(size(E_array) == size(H_array)));

E_array = collapse_complex(E_array);
H_array = collapse_complex(H_array);

% Permute E_array and H_array's indices to (x, y, z, axis)
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
	E_cell{w} = array2scalar(E_array(ind{:}), PhysQ.E, grid3d, w, GK.dual, osc);
	H_cell{w} = array2scalar(H_array(ind{:}), PhysQ.H, grid3d, w, GK.prim, osc);
end
