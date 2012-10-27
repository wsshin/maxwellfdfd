classdef Ellipsoid < Shape
	% Ellipsoid is a Shape for an ellipsoid.

	methods
        function this = Ellipsoid(center, semiaxes, dl_max, dl_boundary)
			chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
				'"center" should be length-%d row vector with real elements.', Axis.count);

			chkarg(istypesizeof(semiaxes, 'real', [1, Axis.count]) && all(semiaxes > 0), ...
				'"semiaxes" should be length-%d row vector with positive elements.', Axis.count);
			
			bound = [center - semiaxes; center + semiaxes];
			bound = bound.';
			
			% The level set function is basically 1 - norm((r-center)./semiaxes)
			% but it is vectorized, i.e., modified to handle r = [x y z] with
			% column vectors x, y, z.
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				N = size(r, 1);
				c_vec = repmat(center, [N 1]);
				s_vec = repmat(semiaxes, [N 1]);
				x = (r - c_vec) ./ s_vec;
				level = 1 - sqrt(sum(x.*x, 2));
			end

			if nargin < 3
				super_args = {bound, @lsf};
			elseif nargin < 4
				super_args = {bound, @lsf, dl_max};
			else
				super_args = {bound, @lsf, dl_max, dl_boundary};
			end
			
			this = this@Shape(super_args);
		end
	end	
end

