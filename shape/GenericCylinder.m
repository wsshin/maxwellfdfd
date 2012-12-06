classdef GenericCylinder < Shape
	% GenericCylinder is a Shape for a cylinder whose cross section is any 2D
	% shape.

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
			sn = diff(bound_n) / 2;
			cn = mean(bound_n);

			% For rhv = r([h, v]), and rn = r(n), the level set
			% function is basically min(lsf2d(rhv), 1 - abs(rn-cn)./sn, [], 2),
			% but it is vectorized, i.e., modified to handle r = [x y z] with
			% column vectors x, y, z.
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				N = size(r, 1);
				cn_vec = cn(ones(N, 1));
				sn_vec = sn(ones(N, 1));
				rhv = r(:, [h v]);
				rn = r(:, n);
				level = min([lsf2d(rhv), 1 - abs(rn-cn_vec)./sn_vec], [], 2);
			end
			
			if nargin < 4  % no dl_max
				super_args = {lprim_cell, @lsf};
			else
				super_args = {lprim_cell, @lsf, dl_max};
			end
			
			this = this@Shape(super_args{:});
		end
	end
end

