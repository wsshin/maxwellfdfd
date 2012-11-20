classdef Line < ZeroVolShape
	% Line is a Shape for a line.
	
	properties
		axis  % direction of line
		intercept  % [h, v]: location of intercept in transverse plane
	end
		
	methods
        function this = Line(axis, intercept, dl_max)
			chkarg(istypesizeof(axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real', [1, Dir.count]), ...
				'"intercept" should be length-%d row vector with real elements.', Dir.count);
			
			bound = [-Inf Inf; -Inf Inf; -Inf Inf];
			[h_axis, v_axis] = cycle(axis);
			bound([h_axis v_axis], :) = [intercept; intercept].';
			
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);

				r_trans = r(:, [h_axis v_axis]);  % transverse position
				N = size(r, 1);
				c = repmat(intercept, [N 1]);
				level = -max(abs(r_trans - c), [], 2);
			end
			
			if nargin < 3  % no dl_max
				super_args = {bound, @lsf};
			else
				dl_max = expand2row(dl_max, Axis.count);
				dl_max(h_axis) = Inf;  dl_max(v_axis) = Inf;  % dl_max is meaningful only in direction of line
				super_args = {bound, @lsf, dl_max};
			end

			this = this@ZeroVolShape(super_args{:});
			this.axis = axis;
			this.intercept = intercept;
		end
	end	
end
