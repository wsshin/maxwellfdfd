classdef DistributedSrc < Source
	% DistributedSrc is a class representing a source distribution on a plane.
	
	properties (SetAccess = immutable)
		normal_axis  % plane normal axis: one of Axis.x, Axis.y, Axis.z
		intercept  % intercept between plane and normal axis
		IL  % (surface current) * (conduction distance); for z-normal plane, volume integral of (|Jx|+|Jy|)
		neff_guess;  % estimated effective refractive index
	end
	
	properties (SetAccess = private)
		grid2d  % instance of Grid2d
		Jh  % J in horizontal direction on this plane: Jx for normal == z
		Jv  % J in vertical direction on this plane: Jy for normal == z
		neff  % effective n
	end
	
	methods
		function this = DistributedSrc(normal_axis, intercept, neff_guess, IL)
			chkarg(istypesizeof(normal_axis, 'Axis'), ...
				'"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			chkarg(istypesizeof(neff_guess, 'complex'), '"neff" should be complex.');
			
			if nargin < 4  % no IL
				IL = 1.0;
			end
			chkarg(istypesizeof(IL, 'real'), '"IL" should be real.');
			
			l = cell(Axis.count, GK.count);
			l{normal_axis, GK.dual} = intercept;
			plane = Plane(normal_axis, intercept);
			this = this@Source(l, plane);
			
			this.normal_axis = normal_axis;
			this.intercept = intercept;
			this.neff_guess = neff_guess;
			this.IL = IL;
		end
		
		function setJ(this, neff, Jh, Jv, grid3d)
			chkarg(istypesizeof(neff, 'complex'), '"neff" should be complex.');
			this.neff = neff;
			
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			this.grid2d = Grid2d(grid3d, this.normal_axis);

			Nh = this.grid2d.N(Dir.h);
			Nv = this.grid2d.N(Dir.v);
			
			assert(istypesizeof(Jh, 'complex', [Nh Nv]), '"Jh" should be %d-by-%d matrix with complex elements.', Nh, Nv);
			assert(istypesizeof(Jv, 'complex', [Nh Nv]), '"Jv" should be %d-by-%d matrix with complex elements.', Nh, Nv);

			this.Jh = Jh;
			this.Jv = Jv;
		end
		
		function [index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
			assert(~isempty(this.Jh) && ~isempty(this.Jv), 'Jh and Jv are not set in this DistributedSrc.');
			index_cell = cell(1, Axis.count);
			if w_axis == this.normal_axis
				Jw_patch = [];
			else
				g2d = Grid2d(grid3d, this.normal_axis);
				assert(isequal(g2d, this.grid2d), ...
					'%s-normal cross section of "grid3d" is different from the one set with Jh and Jv.', char(this.normal_axis));

				h = this.grid2d.axis(Dir.h);
				v = this.grid2d.axis(Dir.v);
				n = this.normal_axis;
				
				g = GK.dual;
				ind_n = ismembc2(this.intercept, grid3d.l{n,g});
				if ind_n == 0
					[~, ind_n] = min(abs(grid3d.l{n,g} - this.intercept));
					warning('FDS:srcAssign', ...
						['%s grid in %s-axis of "grid3d" does not have location %e of this %s; ', ...
						'closest grid vertex at %e will be taken instead.'], ...
						char(g), char(n), this.intercept, class(this), grid3d.l{n,g}(ind_n));
				end
				
				% Set index_cell.
				Nh = this.grid2d.N(Dir.h);
				Nv = this.grid2d.N(Dir.v);
				index_cell{n} = ind_n;
				index_cell{h} = 1:Nh;
				index_cell{v} = 1:Nv;
				
				% Set Jw_patch.
				dn = grid3d.dl{n,g}(ind_n);
				dh_prim = grid3d.dl{h, GK.prim};
				dh_dual = grid3d.dl{h, GK.dual};
				dv_prim = grid3d.dl{v, GK.prim};
				dv_dual = grid3d.dl{v, GK.dual};
				
				dVh = dn .* (dh_prim.' * dv_dual);
				dVv = dn .* (dh_dual.' * dv_prim);
				
				IL_curr = abs(this.Jh) .* dVh + abs(this.Jv) .* dVv;
				IL_curr = sum(IL_curr(:));
				norm_factor = this.IL/IL_curr;  % normalization factor
				
				if w_axis == h
					Jw_patch = norm_factor .* this.Jh;
				else
					assert(w_axis == v);
					Jw_patch = norm_factor .* this.Jv;
				end
				
% 				if h > v  % h == Axis.z and v == Axis.x if this.normal_axis == Axis.y
% 					Jw_patch = permute(Jw_patch, int([Dir.v, Dir.h]));
% 				end
				Jw_patch = ipermute(Jw_patch, int([h v n]));
			end
		end
	end
	
end

