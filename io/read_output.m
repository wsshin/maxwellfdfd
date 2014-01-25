function [E_cell, H_cell] = read_output(filenamebase, israw)
if nargin < 2  % no israw
	israw = false;
end
chkarg(istypesizeof(israw, 'logical'), '"israw" should be logical.');
chkarg(istypesizeof(filenamebase, 'char', [1 0]), '"filenamebase" should be string.');

inputfile = ['./', filenamebase, '.h5'];
efile = ['./', filenamebase, '.E.h5'];
hfile = ['./', filenamebase, '.H.h5'];

lambda = h5read(inputfile, '/lambda');
L0 = h5read(inputfile, '/L0');
Npml = h5read(inputfile, '/Npml').';
ge = h5read(inputfile, '/ge');
ge = GT.elems(ge+1);
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

if israw
	if nD-1 == Axis.count
		E_cell = {E_array(:,:,:,Axis.x), E_array(:,:,:,Axis.y), E_array(:,:,:,Axis.z)};
		H_cell = {H_array(:,:,:,Axis.x), H_array(:,:,:,Axis.y), H_array(:,:,:,Axis.z)};
	else
		assert(nD-1 == Dir.count);
		E_cell = {E_array(:,:,Axis.x), E_array(:,:,Axis.y), E_array(:,:,Axis.z)};
		H_cell = {H_array(:,:,Axis.x), H_array(:,:,Axis.y), H_array(:,:,Axis.z)};
	end
else
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
end
