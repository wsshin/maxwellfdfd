function [Dive, Divm] = create_divs(ge, dl_cell, grid3d)

chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for the E-field
chkarg(istypesizeof(dl_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
	'"dl_cell" should be %d-by-%d cell array whose each element is row vector with real elements.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

%% Create Df and Db.
Df = create_Ds(Sign.p, ge, dl_cell, grid3d);
Db = create_Ds(Sign.n, ge, dl_cell, grid3d);

%% Create mask matrices.
% [Mp, Md] = create_masks(ge, grid3d);

%% Form matrices
Divp = [Db{Axis.x}, Db{Axis.y}, Db{Axis.z}];
% Divp = Divp * Mp;

Divd = [Df{Axis.x}, Df{Axis.y}, Df{Axis.z}];
% Divd = Divd * Md;

if ge == GT.prim
	Dive = Divp;
	Divm = Divd;
else  % ge == GT.dual
	Dive = Divd;
	Divm = Divp;
end
