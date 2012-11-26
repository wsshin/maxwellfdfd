classdef Rectangle < ZeroVolShape
	% Plane is a Shape for a plane.
	
	properties
		normal_axis
		intercept
		rect
	end
		
	methods
        function this = Rectangle(normal_axis, intercept, rect, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			chkarg(istypesizeof(rect, 'real', [Dir.count Sign.count]), ...
				'"rect" should be %d-by-%d array with real elements.', Dir.count, Sign.count);

			bound = NaN(Axis.count, Sign.count);
			bound(normal_axis, :) = [intercept intercept];
			[h, v] = cycle(normal_axis);
			bound([h v], :) = rect;
			
			box = Box(bound);
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				r_normal = r(:, normal_axis);
				r(:, normal_axis) = intercept;  % without this, box.lsf() returns -Inf in most cases, where ZeroVolShape.lsf(r, true) is in vain
				level = box.lsf(r) - abs(r_normal - intercept);
			end
			
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				lprim{w} = box.bound(w,:);
			end
					
			if nargin < 4  % no dl_max
				super_args = {lprim, @lsf};
			else
				dl_max = expand2row(dl_max, Axis.count);
				dl_max(normal_axis) = Inf;  % dl_max is meaningless in normal direction
				super_args = {lprim, @lsf, dl_max};
			end

			this = this@ZeroVolShape(super_args{:});
			this.normal_axis = normal_axis;
			this.intercept = intercept;
			this.rect = rect;
		end
	end	
end
