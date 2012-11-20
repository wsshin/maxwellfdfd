classdef PointSrc < Source
	% PointSrc is a class representing an electric point source.
	
	properties (SetAccess = immutable)
		polarization  % one of Axis.x, Axis.y, Axis.z
		location  % [x, y, z]
		I  % current, e.g., Jz * dx * dy
	end
	
	methods
		function this = PointSrc(polarization_axis, location, I)
			chkarg(istypesizeof(polarization_axis, 'Axis'), ...
				'"polarization_axis" should be instance of Axis.');
			chkarg(istypesizeof(location, 'real', [1, Axis.count]), ...
				'"location" should be length-%d row vector with real elements.', Axis.count);
			if nargin < 3  % no I
				I = 1.0;
			end
			chkarg(istypesizeof(I, 'complex'), '"I" should be complex.');
			
			l = cell(Axis.count, GK.count);
			for w = Axis.elems
				if w == polarization_axis
					l{w, GK.prim} = location(w);
				else
					l{w, GK.dual} = location(w);
				end
			end
			point = Point(location);
			this = this@Source(l, point);
			this.polarization = polarization_axis;
			this.location = location;
			this.I = I;
		end
				
		function [index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
			index_cell = cell(1, Axis.count);
			[q, r, p] = cycle(this.polarization);  % q, r: axes normal to polarization axis p
			if w_axis == p
				for v = Axis.elems
					l = this.location(v);
					if v == p
						g = GK.prim;
					else  % v == q or r
						g = GK.dual;
					end
					iv = ismembc2(l, grid3d.l{v,g});
					if iv == 0
						[~, iv] = min(abs(grid3d.l{v,g} - l));
						warning('FDS:srcAssign', ...
							['%s grid in %s-axis of "grid3d" does not have location %e of this %s; ', ...
							'closest grid vertex at %e will be taken instead.'], ...
							char(g), char(v), l, class(this), grid3d.l{v,g}(iv));
					end
					index_cell{v} = iv;
				end
				dq = grid3d.dl{q, GK.dual}(index_cell{q});
				dr = grid3d.dl{r, GK.dual}(index_cell{r});
				Jw_patch = this.I / (dq * dr);  % I = J * (area)
			else  % w_axis == q or r
				Jw_patch = [];
			end
		end		
	end
end
