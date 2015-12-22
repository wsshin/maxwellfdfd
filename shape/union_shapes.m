function lsf = union_shapes(ind2shape_array, shape_ind_array)
	chkarg(istypesizeof(ind2shape_array, 'Shape', [1 0]), '"ind2shape_array" is row vector with Shape as elements.');
	chkarg(istypesizeof(shape_ind_array, 'int', [1 0]), '"shape_ind_array" is row vector with integer elements.');
	
	function level = lsf_union(x, y, z)
		chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
		chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
		chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
		chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');
		
		level = -Inf(size(x));
		for s_ind = shape_ind_array
			s = ind2shape_array(s_ind);
			level = max(level, s.lsf(x, y, z));
		end
	end

	lprim = cell(1, Axis.count);  % empty
	for shape_ind = shape_ind_array
		shape = ind2shape_array(shape_ind);
		for w = Axis.elems
			lprim{w} = [lprim{w}, shape.lprim{w}];
		end
	end
	lsf = Shape(lprim, @lsf_union);	
end
