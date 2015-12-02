classdef GenericCylinder < Shape
	% GenericCylinder is a Shape for a cylinder whose cross section is any 2D
	% shape.
	
	properties (SetAccess = immutable)
		lsf2d  % level set function in lateral domain
	end

	methods
        function this = GenericCylinder(normal_axis, lsf2d, lprim_cell, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			
			% lsf2d is a level set function defined in 2D.  It takes an argument
			% r = [h, v], where, e.g., (h, v) is (x, y) for narmal_axis: z.
			chkarg(istypeof(lsf2d, 'function_handle'), '"lsf2d" should be function handle.');

			chkarg(istypesizeof(lprim_cell, 'realcell', [1 Axis.count], [1 0]), ...
				'"lprim_cell" should be length-%d row cell array whose each element is row vector with real elements.', Axis.count);
			
			[h, v, n] = cycle(normal_axis);
			bound_n = [min(lprim_cell{n}), max(lprim_cell{n})];
			sn = diff(bound_n) / 2;  % semiside in normal direction
			cn = mean(bound_n);  % center in normal direction

			% For rhv = r([h, v]), and rn = r(n), the level set
			% function is basically min(lsf2d(rhv), 1 - abs(rn-cn)./sn, [], 2),
			% but it is vectorized, i.e., modified to handle r = [x y z] with
			% column vectors x, y, z.
			function level = lsf(x, y, z)
				chkarg(istypeof(x, 'real'), '"x" should be array with real elements.');
				chkarg(istypeof(y, 'real'), '"y" should be array with real elements.');
				chkarg(istypeof(z, 'real'), '"z" should be array with real elements.');
				chkarg(isequal(size(x), size(y), size(z)), '"x", "y", "z" should have same size.');
				
				loc = {x, y, z};
				level = min(lsf2d(loc{h}, loc{v}), 1 - abs(loc{n}-cn)./sn);  % intersection of regions defined by two level set functions
			end
			
			if nargin < 4  % no dl_max
				super_args = {lprim_cell, @lsf};
			else
				super_args = {lprim_cell, @lsf, dl_max};
			end
			
			this = this@Shape(super_args{:});
			this.lsf2d = lsf2d;
		end
	end
end

