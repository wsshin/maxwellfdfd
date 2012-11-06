classdef PointSrcM < Source
	% PointSrc is a class representing an magnetic point source.
	
	properties (SetAccess = immutable)
		polarization  % one of Axis.x, Axis.y, Axis.z
		location  % [x, y, z];
		Im  % current, e.g., Jz * dx * dy
	end
	
	methods
		function this = PointSrcM(polarization_axis, location, Im)
			chkarg(istypesizeof(polarization_axis, 'Axis'), ...
				'"polarization_axis" should be instance of Axis.');
			chkarg(istypesizeof(location, 'real', [1, Axis.count]), ...
				'"location" should be length-%d row vector with real elements.', Axis.count);
			if nargin < 3  % no I
				Im = 1.0;
			end
			chkarg(istypesizeof(Im, 'complex'), '"I" should be complex.');
			
			l = cell(Axis.count, GK.count);
			for w = Axis.elems
				if w == polarization_axis
					l{w, GK.dual} = location(w);
				else
					l{w, GK.prim} = location(w);
				end
			end
			this = this@Source(l);
			this.polarization = polarization_axis;
			this.location = location;
			this.Im = Im;
		end
				
		function [index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
			index_cell = cell(1, Axis.count);
			[q, r, p] = cycle(this.polarization);  % q, r: axes normal to polarization axis p
			if w_axis == p
				Jw_patch = [];
			else  % w_axis == q or r
				indM = NaN(1, Axis.count);
				for v = Axis.elems
					l = this.location(v);
					if v == p
						g = GK.dual;
					else  % v == q or r
						g = GK.prim;
					end
					iv = ismembc2(l, grid3d.l{v,g});
					if iv == 0
						[~, iv] = min(abs(grid3d.l{v,g} - l));
						warning('FDS:srcAssign', ...
							['%s grid in %s-axis of "grid3d" does not have location %e of this %s; ', ...
							'closest grid vertex at %e will be taken instead.'], ...
							char(g), char(v), l, class(this), grid3d.l{v,g}(iv));
					end
					indM(v) = iv;
				end
				dq = grid3d.dl{q, GK.prim}(indM(q));
				dr = grid3d.dl{r, GK.prim}(indM(r));
				I = this.Im / (dq * dr);  % (magnetic dipole) = (electric current) * (area)
				
				
				if w_axis == q
					% Assign Jq.
					index_cell{p} = indM(p);
					index_cell{q} = indM(q);
					index_cell{r} = [indM(r) - 1, indM(r)];
					assert(indM(r) - 1 >= 1, ...
						'PointSrcM should not be at boundary of %s-axis.', char(r));
					dlp = grid3d.dl{r, GK.dual}(index_cell{p});  % one element
					dlr = grid3d.dl{r, GK.dual}(index_cell{r});  % two elements
					dS = dlp .* dlr;
					Jw_patch = [I, -I] ./ dS;
				else
					assert(w_axis == r);
					% Assign Jr.
					index_cell{p} = indM(p);
					index_cell{q} = [indM(q) - 1, indM(q)];
					index_cell{r} = indM(r);
					assert(indM(q) - 1 >= 1, ...
						'PointSrcM should not be at boundary of %s-axis.', char(q));
					dlp = grid3d.dl{r, GK.dual}(index_cell{p});  % one element
					dlq = grid3d.dl{q, GK.dual}(index_cell{q});  % two elements
					dS = dlp .* dlq;
					Jw_patch = [-I, I] ./ dS;
					Jw_patch = Jw_patch.';
				end
				Jw_patch = ipermute(Jw_patch, int([q r p]));
			end
		end		
	end
end
