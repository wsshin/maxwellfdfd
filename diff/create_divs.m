function [Dive, Divm] = create_divs(ge, dl_factor_cell, grid3d)

chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for the E-field
chkarg(isempty(dl_factor_cell) || istypesizeof(dl_factor_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
	'"dl_factor_cell" should be empty, or %d-by-%d cell array whose each element is row vector with real elements.', Axis.count, GT.count);
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

%% Create Df and Db.
% Surprisingly, it seems that Ds used to create the curl operators can be reused
% to create the divergence operators, while satisfying all the boundary
% conditions.
Df = create_Ds(Sign.p, ge, dl_factor_cell, grid3d);
Db = create_Ds(Sign.n, ge, dl_factor_cell, grid3d);

%% Create mask matrices.
% [Mp, Md] = create_masks(ge, grid3d);

%% Form matrices
Divp = [Db{Axis.x}, Db{Axis.y}, Db{Axis.z}];  % primary fields use backward difference
% Divp = Divp * Mp;

Divd = [Df{Axis.x}, Df{Axis.y}, Df{Axis.z}];  % dual fields use forward difference
% Divd = Divd * Md;

if ge == GT.prim
	Dive = Divp;
	Divm = Divd;
else  % ge == GT.dual
	Dive = Divd;
	Divm = Divp;
end
