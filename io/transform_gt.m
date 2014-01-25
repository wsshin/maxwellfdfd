function scalar = transform_gt(scalarNd, axis, gt_array)

chkarg(istypesizeof(scalarNd, 'Scalar2d') || istypesizeof(scalarNd, 'Scalar3d'), ...
	'"scalarNd" should be instance of Scalar2d or Scalar3d.');
chkarg(istypesizeof(axis, 'Axis'), '"axis" should be instance of Axis.');

physQcell = scalarNd.physQcell;
chkarg(length(physQcell)==2 && physQcell{2}==1, '"scalarNd" should represent quantity defined in PhysQ.');

physQ = scalarNd.physQcell{1};
chkarg(physQ==PhysQ.E || physQ==PhysQ.H, '"physQ" should be either PhysQ.E or PhysQ.H.');
if physQ == PhysQ.E
	ft = FT.e;
else
	ft = FT.h;
end

if istypesizeof(scalarNd, 'Scalar2d')
	Ndim = 2;
	grid = scalarNd.grid2d;
	axis_dummy = Dir.h;
else
	Ndim = 3;
	grid = scalarNd.grid3d;
	axis_dummy = Axis.x;
end

chkarg(istypesizeof(gt_array, 'GT', [1 Ndim]), ...
	'"gt_array" should be length-%d row vector of GT objects.', Ndim);

[data, l] = scalarNd.data_expanded();
li = grid.l(axis_dummy.elems + axis_dummy.count*subsindex(gt_array));  % locations where data are interpreted

L = cell(1,Ndim);
LI = cell(1,Ndim);
[L{:}] = ndgrid(l{:});
[LI{:}] = ndgrid(li{:});

input = horzcat(L, {data}, LI);
V = interpn(input{:});

osc = scalarNd.osc;

if Ndim == 3
	scalar = array2scalar(V, physQ, grid, axis, ft, gt_array, osc);
else  % Ndim ==2
	intercept = scalarNd.intercept;
	scalar = array2scalar(V, physQ, grid, axis, ft, gt_array, osc, intercept);
end
