function patch_handle_array = draw_objsrc(obj_array, grid3d, withinterp, withpml)

chkarg(istypesizeof(obj_array, 'Object', [1 0]), ...
	'"object_array" should be row vector with Object as elements.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');

if withinterp
	lplot = grid3d.lplot(GT.prim, withinterp, withpml);
else
	lplot = grid3d.lpixelbound(GT.prim, withpml);
end

for w = Axis.elems
	if length(lplot{w}) <= 3
		lplot{w} = linspace(lplot{w}(1), lplot{w}(end), 10);
	end
end

boxes_pml = Object.empty(0,0);
gray = [0.5 0.5 0.5];
if withpml
	for w = Axis.elems
		for s = Sign.elems
			if grid3d.Npml(w,s) > 0
				bound = grid3d.bound;
				bound(w, alter(s)) = grid3d.lpml(w,s);
				box = Box(bound);
				pml = Material('PML', gray, 1.0);
				box_pml = Object(box, pml);
				boxes_pml = [boxes_pml(1:end), box_pml];
			end
		end
	end
end

patch_handle_array = [];
lplotobj = cell(1, Axis.count);  % indices
for obj = [obj_array, boxes_pml]
	shape = obj.shape;
	color = obj.material.color;
	if ~isequal(color, 'none') && ~istypesizeof(shape, 'Domain')
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
			lsf = shape.lsf;
			[X, Y, Z] = meshgrid(lplotobj{:});
			Level = lsf([X(:), Y(:), Z(:)]);
			Level = reshape(Level, size(X));
			hp = patch(isosurface(X, Y, Z, Level, 0));
			isonormals(X, Y, Z, Level, hp)

			set(hp, 'FaceColor', color, 'EdgeColor', 'none');
			if isequal(obj.material.name, 'PML')
				alpha(hp, 0.5);
			end

			patch_handle_array = [patch_handle_array(1:end), hp];
		end
	end
end


