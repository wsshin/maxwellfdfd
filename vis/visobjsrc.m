function visobjsrc(grid3d, obj_array, src_array, withpml)

chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(obj_array, 'Object', [1 0]), '"obj_array" should be row vector with Object as elements.');
chkarg(istypesizeof(src_array, 'Source', [1 0]), '"src_array" should be row vector with Source as elements.');

if nargin < 3  % no withpml
	withpml = true;
end
chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');

init_axes3d(gca, grid3d, true, withpml);

% Below, pass "withinterp = false" for 2D simulation.  Because no primary grid
% point exists inside a shape in 3rd direction in 2D simulation, with
% "withinterp = true", which examines primary grid nodes, no object is
% visualized.
draw_objsrc(obj_array, grid3d, false, withpml);  
camlight
lighting gouraud
