function [Ce, Cm] = create_curls(ge, dl_factor_cell, grid3d)

chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for the E-field
chkarg(isempty(dl_factor_cell) || istypesizeof(dl_factor_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
	'"dl_factor_cell" should be empty, or %d-by-%d cell array whose each element is row vector with real elements.', Axis.count, GT.count);
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');


%% Create Df and Db.
Df = create_Ds(Sign.p, ge, dl_factor_cell, grid3d);
Db = create_Ds(Sign.n, ge, dl_factor_cell, grid3d);

%% Create mask matrices.
[ind_Mp, ind_Md] = create_masks(ge, grid3d);

%% Form curl matrices.
M = prod(grid3d.N);
Z = sparse(M, M);

Cp = [Z, -Df{Axis.z}, Df{Axis.y}; ...
	Df{Axis.z}, Z, -Df{Axis.x}; ...
	-Df{Axis.y}, Df{Axis.x}, Z];
% Cp = Md * Cp * Mp;
Cp(:, ind_Mp) = 0;
Cp(ind_Md, :) = 0;

Cd = [Z, -Db{Axis.z}, Db{Axis.y}; ...
	Db{Axis.z}, Z, -Db{Axis.x}; ...
	-Db{Axis.y}, Db{Axis.x}, Z];
% Cd = Mp * Cd * Md;
Cd(:, ind_Md) = 0;
Cd(ind_Mp, :) = 0;

if ge == GT.prim
	Ce = Cp;
	Cm = Cd;
else  % ge == GT.dual
	Ce = Cd;
	Cm = Cp;
end
