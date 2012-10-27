function patch_handle_array = draw_objsrc(objec_array, grid3d, withinterp, withpml)

chkarg(istypesizeof(objec_array, 'Object', [1 0]), ...
	'"object_array" should be row vector with Object as elements.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');

patch_handle_array = [];
if withinterp
	lplot = grid3d.lplot(GK.prim, withpml);
else
	lplot = grid3d.lplot(GK.dual, withpml);
end

ind = cell(1, Axis.count);  % indices
for obj = objec_array
	color = obj.material.color;
	if ~isequal(color, 'none')
		shape = obj.shape;
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
			ind{w} = in:ip;
		end

		if ~isempty(ind{Axis.x}) && ~isempty(ind{Axis.y}) && ~isempty(ind{Axis.z})  % shape in inside domain
			lsf = shape.lsf;
			[X, Y, Z] = meshgrid(lplot{Axis.x}(ind{Axis.x}), lplot{Axis.y}(ind{Axis.y}), lplot{Axis.z}(ind{Axis.z}));
			Level = lsf([X(:), Y(:), Z(:)]);
			Level = reshape(Level, size(X));
			hp = patch(isosurface(X, Y, Z, Level, 0));
			isonormals(X, Y, Z, Level, hp)

			set(hp, 'FaceColor', color, 'EdgeColor', 'none');

			patch_handle_array = [patch_handle_array(1:end), hp];
		end
	end
end

