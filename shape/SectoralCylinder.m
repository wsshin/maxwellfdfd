classdef SectoralCylinder < GenericCylinder
	% EllipticCylinder is a Shape for a cylinder whose cross section is an
	% ellipse.

	methods
        function this = SectoralCylinder(normal_axis, height, center, radius, theta, d_theta, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');
			chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
				'"center" should be length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(radius, 'real') && radius > 0, '"radius" should be positive.');
			chkarg(istypesizeof(theta, 'real'), '"theta" should be real.');
			chkarg(istypesizeof(d_theta, 'real') && (d_theta <= 2*pi || d_theta >= -2*pi), '"d_theta" should be real between -pi and pi.');
			
			thetas = sort([theta, theta + d_theta]);
			cth = mean(thetas);
			cth = mod(cth + pi, 2*pi) - pi;  % -pi <= c_th < pi (range of atan2(y,x))
			sth = diff(thetas) / 2;
			
			cr = radius/2;
			sr = radius/2;
			
			function level = lsf_th(th)
				dth = th - cth;
				dth = mod(dth + pi, 2*pi) - pi;
				level = 1 - abs(dth/sth);
			end

			[h, v, n] = cycle(normal_axis);

			% lsf2d() can handle rho = [p q] with column vectors p and q.  The
			% level set function is the one for a rectangle defined in the
			% (theta, radius) domain.
			function level = lsf2d(rho)
				chkarg(istypesizeof(rho, 'real', [0, Dir.count]), ...
					'"rho" should be matrix with %d columns with real elements.', Dir.count);
				N = size(rho, 1);
				c = center([h, v]);
				c_vec = repmat(c, [N 1]);
				
				rc = rho - c_vec;
				theta_pt = atan2(rc(:,Dir.v), rc(:,Dir.h));  % -pi <= theta_rho < pi
				r_pt = sqrt(sum(rc.^2, 2));
				zero_r_pt = (r_pt==0);  % atan2(0,0) is not well-defined
				
				dr = r_pt - cr;
				level = min([lsf_th(theta_pt), 1 - abs(dr/sr)], [], 2);
				level(zero_r_pt) = 0;
			end
			
			lprim = cell(1, Axis.count);
			lprim{h} = [radius * cos(thetas) + center(h), center(h)];
			lprim{v} = [radius * sin(thetas) + center(v), center(v)];
			lprim{n} = [-height/2 height/2] + center(n);
			
			if lsf_th(0) > 0
				lprim{h} = [lprim{h}, center(h) + radius];
			end
			if lsf_th(pi/2) > 0
				lprim{v} = [lprim{v}, center(v) + radius];
			end
			if lsf_th(pi) > 0
				lprim{h} = [lprim{h}, center(h) - radius];
			end
			if lsf_th(3*pi/2) > 0
				lprim{v} = [lprim{v}, center(v) - radius];
			end
						
			if nargin < 7  % no dl_max
				super_args = {normal_axis, @lsf2d, lprim};
			else
				super_args = {normal_axis, @lsf2d, lprim, dl_max};
			end
			
			this = this@GenericCylinder(super_args{:});
		end
	end
end

