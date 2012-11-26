classdef Point < ZeroVolShape
	% Line is a Shape for a line.
	
	properties
		location  % [x, y, z]: location of point
	end
		
	methods
        function this = Point(location, dl_max)
			chkarg(istypesizeof(location, 'real', [1, Axis.count]), ...
				'"intercept" should be length-%d row vector with real elements.', Axis.count);
			
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				
				N = size(r, 1);
				c = repmat(location, [N 1]);
				level = -max(abs(r - c), [], 2);
			end
			
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				lprim{w} = location(w);
			end
			
			super_args = {lprim, @lsf};  % dl_max is meaningless

			this = this@ZeroVolShape(super_args{:});
			this.location = location;
		end
	end	
end
