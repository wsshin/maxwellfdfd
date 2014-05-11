function [grid3d, ge, osc] = read_grid(filenamebase)

chkarg(istypesizeof(filenamebase, 'char', [1 0]), '"filenamebase" should be string.');

inputfile = ['./', filenamebase, '.h5'];

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
