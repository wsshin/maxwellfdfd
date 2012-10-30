classdef Box < Shape
	% Box is a Shape for a box.
		
	methods
        function this = Box(bound, dl_max)
			% bound
			chkarg(istypesizeof(bound, 'real', [Axis.count, Sign.count]), ...
				'"bound" should be [xmin xmax; ymin ymax; zmin zmax].');
			s = diff(bound, 1, 2) ./ 2;  % semisides
			s = s.';  % row vector
			chkarg(all(s > 0), '"bound" should have smaller lower bound than upper bound in all axes.');
			c = mean(bound, 2);  % center
			c = c.';  % row vector
			
			% The level set function is basically 1 - max(abs(r-c) ./ s, [], 2),
			% but it is vectorized, i.e., modified to handle r = [x y z] with
			% column vectors x, y, z.
			function level = lsf(r)
				chkarg(istypesizeof(r, 'real', [0, Axis.count]), ...
					'"r" should be matrix with %d columns with real elements.', Axis.count);
				N = size(r, 1);
				c_vec = repmat(c, [N 1]);
				s_vec = repmat(s, [N 1]);
				level = 1 - max(abs(r - c_vec) ./ s_vec, [], 2);
			end
			
			if nargin < 2  % no dl_max
				super_args = {bound, @lsf};
			else
				super_args = {bound, @lsf, dl_max};
			end

			this = this@Shape(super_args{:});
		end
		
% 		function [n_dir, r_vol] = ndir_and_rvol(this, pixel)
% 			chkarg(istypesizeof(pixel, 'real', [Axis.count, Sign.count]), '"pixel" should be [xmin xmax; ymin ymax; zmin zmax].');
% 			center = mean(pixel, 2).';
% 			
% 			in = NaN(1, Axis.count);
% 			d = NaN(1, Axis.count);
% 			for w = Axis.elems
% 				[in(w), d(w)] = this.interval(w).contains(center(w));
% 			end
% 			chkarg((in(Axis.x) && in(Axis.y)) || (in(Axis.y) && in(Axis.z)) || (in(Axis.z) && in(Axis.x)), 'center of "pixel" should be contained the box at least in two directions.');
% 
% 			if all(in)
% 				[~, i] = min(d);
% 			else  % only one element of "in" is false
% 				i = find(~in);
% 			end
% 			
% 			w = Axis.elems(i);
% 			n_dir = zeros(1,Axis.count);
% 			n_dir(w) = 1;
% 			dc = d(w);
% 			L = diff(pixel(w,:));
% 			chkarg(L, '"pixel" should be [xmin xmax; ymin ymax; zmin zmax], where wmin < wmax.');
% 			if dc >= L/2
% 				r_vol = 0;
% 			elseif in(w)
% 				r_vol = (L/2 + dc) / L;
% 			else  % in(w) == false
% 				r_vol = (L/2 - dc) / L;
% 			end
% 		end
	end	
end
