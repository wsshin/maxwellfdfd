classdef UnionShape < Shape

	properties (SetAccess = immutable)
		shape_array  % array of component shapes
	end
	
	methods
        function this = UnionShape(shape_array)
			chkarg(istypesizeof(shape_array, 'Shape', [1 0]), '"shape_array" is row vector with Shape as elements.');

			function level = lsf_union(x, y, z)
				chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
				chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
				chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
				chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');

				level = -Inf(size(x));
				for s = shape_array
					level = max(level, s.lsf(x, y, z));
				end
			end

			lprim = cell(1, Axis.count);  % empty
			for shape = shape_array
				for w = Axis.elems
					lprim{w} = [lprim{w}, shape.lprim{w}];
				end
			end
			
			super_args = {lprim, @lsf_union};
			this = this@Shape(super_args{:});
			
			this.shape_array = shape_array;
		end
	end
end
