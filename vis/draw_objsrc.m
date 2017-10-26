function patch_handle_array = draw_objsrc(obj_array, src_array, grid3d, withinterp, withpml)

chkarg(istypesizeof(obj_array, 'EMObject', [1 0]), ...
	'"obj_array" should be row vector with EMObject as elements.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');

if withinterp
	lplot = grid3d.lplot(GT.prim, withinterp, withpml);
else
	lplot = grid3d.lvoxelbound(GT.prim, withpml);
end

for w = Axis.elems
	if length(lplot{w}) <= 3
		lplot{w} = linspace(lplot{w}(1), lplot{w}(end), 10);
	end
end

boxes_pml = EMObject.empty(0,0);
gray = [0.5 0.5 0.5];
if withpml
	for w = Axis.elems
		for s = Sign.elems
			if grid3d.Npml(w,s) > 0
				bound = grid3d.bound;
				bound(w, alter(s)) = grid3d.lpml(w,s);
				box = Box(bound);
				pml = Material('PML', gray, 1.0);
				box_pml = EMObject(box, pml);
				boxes_pml = [boxes_pml(1:end), box_pml];
			end
		end
	end
end

srcobj_array = EMObject.empty(0, length(src_array));
i = 0;
for src = src_array
	i = i+1;
	green = 'g';  % green
	srcmat = Material('Source', green, 1.0);
	srcobj_array(i) = EMObject(src.shape, srcmat);
end

patch_handle_array = [];
lplotobj = cell(1, Axis.count);  % indices
for obj = [obj_array, srcobj_array, boxes_pml]
	shape = obj.shape;
	color = obj.material.color;
% 	if ~isequal(color, 'none') && ~istypesizeof(shape, 'Domain')
	if ~isequal(color, 'none')
		for w = Axis.elems
			bn = shape.bound(w, Sign.n);
			bp = shape.bound(w, Sign.p);
			in = find(lplot{w} >= bn, 1, 'first');
			ip = find(lplot{w} <= bp, 1, 'last');
			
			% Find the indices for a box that is slighly larger than "shape".
			if in > 1
				in = in - 1;
			end
			if ip < length(lplot{w})
				ip = ip + 1;
			end
			lw = lplot{w}(in:ip);
			Nmin = 20;  % minimun number of sampling points
			if ~isempty(lw) && length(lw) < Nmin
				lw = linspace(lw(1), lw(end), Nmin);
			end
			lplotobj{w} = lw;
		end

		if ~isempty(lplotobj{Axis.x}) && ~isempty(lplotobj{Axis.y}) && ~isempty(lplotobj{Axis.z})  % shape in inside domain
			if ~istypesizeof(shape, 'ZeroVolShape')
				lsf = shape.lsf;
			else
				lsf = @(x,y,z) shape.lsf(x, y, z, true);
			end
			[X, Y, Z] = meshgrid(lplotobj{:});
			Level = lsf(X, Y, Z);
			hp = patch(isosurface(X, Y, Z, Level, -eps));  % -eps instead of 0 for ZeroVolumeShape
			isonormals(X, Y, Z, Level, hp)

			if istypesizeof(shape, 'Point')
				set(hp, 'FaceColor', color, 'EdgeColor', color, 'LineWidth', 3);
			elseif istypesizeof(shape, 'Point') || istypesizeof(shape, 'Plane')
				set(hp, 'FaceColor', color, 'EdgeColor', color, 'LineWidth', 1);
			else
				set(hp, 'FaceColor', color, 'EdgeColor', 'none');
			end
			if isequal(obj.material.name, 'PML') || isequal(obj.material.name, 'Source')
				alpha(hp, 0.2);
			end

			patch_handle_array = [patch_handle_array(1:end), hp];
		end
	end
end


