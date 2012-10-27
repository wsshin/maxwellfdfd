classdef EllipticCylinder < GenericCylinder
	% EllipticCylinder is a Shape for a cylinder whose cross section is an
	% ellipse.

	methods
        function this = EllipticCylinder(normal_axis, center, semiaxes, height, dl_max, dl_boundary)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
				'"center" should be length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(semiaxes, 'real', [1, Dir.count]) && all(semiaxes > 0), ...
				'"semiaxes" should be length-%d row vector with positive elements.', Axis.count);
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');

			[h v n] = cycle(normal_axis);
			s = NaN(1, Axis.count);  % semisides
			s(h) = semiaxes(Dir.h);
			s(v) = semiaxes(Dir.v);
			s(n) = height / 2;
			bound = [center - s; center + s];
			bound = bound.';
			
			% For c = center([h, v]), the level set function is
			% basically 1 - norm((rho-c) ./ semiaxes) but it is vectorized,
			% i.e., modified to handle rho = [p q] with column vectors p and q.
			function level = lsf2d(rho)
				chkarg(istypesizeof(rho, 'real', [0, Dir.count]), ...
					'"rho" should be matrix with %d columns with real elements.', Dir.count);
				N = size(rho, 1);
				c = center([h, v]);
				c_vec = repmat(c, [N 1]);
				s_vec = repmat(semiaxes, [N 1]);
				x = (rho - c_vec) ./ s_vec;
				level = 1 - sqrt(sum(x.*x, 2));
			end
			
			if nargin < 5  % no dl_max
				super_args = {normal_axis, @lsf2d, bound};
			elseif nargin < 6  % no dl_boundary
				super_args = {normal_axis, @lsf2d, bound, dl_max};
			else
				super_args = {normal_axis, @lsf2d, bound, dl_max, dl_boundary};
			end
			
			this = this@GenericCylinder(super_args{:});
		end
	end
end

