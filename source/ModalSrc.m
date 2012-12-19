%% ModalSrc
% Concrete subclass of <Source.html |Source|> representing an electric dipole
% distribution that generates a mode of a waveguide.

%%% Description
% |ModalSrc| is used to store the transverse _J_-field distribution on a plane
% to generate a mode of a waveguide.  The plane should be orthogonal to the axis
% of the waveguide.  The mode is calculated by an external mode solver.  To
% assist the mode solver to calculate the mode, users need to provide an
% estimate of the effective refractive constant of the mode to the constructor
% of |ModalSrc|.

%%% Construction
%  src = ModalSrc(normal_axis, intercept, neff_guess)
%  src = ModalSrc(normal_axis, intercept, neff_guess, KA)
% 
% *Input Arguments*
%
% * |normal_axis|: direction normal to the plane over which the electric dipoles
% distribute.  It is also the axis of the waveguide.  It should be one of
% |Axis.x|, |Axis.y|, |Axis.z|.
% * |intercept|: location of the plane in the |normal_axis| direction.
% * |neff_guess|: estimate of the effective refractive index of the mode to
% calculate. 
% * |KA|: integral of the norm of the surface current density.  It indicates the
% strength of the mode amplitude.  If not assigned, the deault value |KA = 1| is
% used.

%%% Note
% In the finite-difference grid, |ModalSrc| excites dipoles at the _E_-field
% points.  This poses a condition on |intercept| argument in the constructor:
% |intercept| should be at a dual grid point in the |normal_axis| direction.
% Therefore, make sure that |intercept| does not overlap with the locations of
% the vertices of <Shape.html |Shape|> in the |normal_axis| direction; otherwise
% dynamic grid generation in <moxwell_run.html |maxwell_run|> will fail.

%%% Example
%   % Create an instance of PointSrc.
%   src =  ModalSrc(Axis.y, -1000, 1.0);  % y = -1000 should not be primary grid point
%
%   % Use the constructed src in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'SRC', src);

%%% See Also
% <Plane.html |PlaneSrc|>, <TFSFPlaneSrc.html |TFSFPlaneSrc|>, <ModalSrc.html
% |ModalSrc|>, <maxwell_run.html |maxwell_run|>

classdef ModalSrc < Source
	
	properties (SetAccess = immutable)
		normal_axis  % plane normal axis: one of Axis.x, Axis.y, Axis.z
		intercept  % intercept between plane and normal axis
		KA  % surface integral of surface current density K; for z-normal plane, volume integral of (|Jx|+|Jy|)
		neff_guess;  % estimated effective refractive index
	end
	
	properties (SetAccess = private)
		grid2d  % instance of Grid2d
		Ft2d  % {Hh2d Hv2d Hn2d; Eh2d Ev2d En2d}: cell array of Scalar2d for H and E on this plane; {Hx2d Hy2d Hz2d; Ex2d Ey2d Ez2d} for normal == z
		Jh  % 2D array J in horizontal direction on this plane: Jx for normal == z
		Jv  % 2D array J in vertical direction on this plane: Jy for normal == z
		neff  % effective n
	end
	
	methods
		function this = ModalSrc(normal_axis, intercept, neff_guess, KA)
			chkarg(istypesizeof(normal_axis, 'Axis'), ...
				'"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			chkarg(istypesizeof(neff_guess, 'complex'), '"neff" should be complex.');
			
			if nargin < 4  % no KA
				KA = 1.0;
			end
			chkarg(istypesizeof(KA, 'real'), '"KA" should be real.');
			
			l = cell(Axis.count, GK.count);
			l{normal_axis, GK.dual} = intercept;
			plane = Plane(normal_axis, intercept);
			this = this@Source(l, plane);
			
			this.normal_axis = normal_axis;
			this.intercept = intercept;
			this.neff_guess = neff_guess;
			this.KA = KA;
			this.Jh = [];
			this.Jv = [];
		end
		
		function setEH(this, neff, osc, E_cell, H_cell, grid3d)
			chkarg(istypesizeof(neff, 'complex'), '"neff" should be complex.');
			this.neff = neff;
			
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			this.grid2d = Grid2d(grid3d, this.normal_axis);

			Nh = this.grid2d.N(Dir.h);
			Nv = this.grid2d.N(Dir.v);
			h = this.grid2d.axis(Dir.h);
			v = this.grid2d.axis(Dir.v);
			
			assert(istypesizeof(E_cell, 'complexcell', [1 Axis.count], [Nh Nv]), ...
				'"E_cell" should be length-%d row cell array whose each element is %d-by-%d array with complex elements.', Axis.count, Nh, Nv);
			assert(istypesizeof(H_cell, 'complexcell', [1 Axis.count], [Nh Nv]), ...
				'"H_cell" should be length-%d row cell array whole each element is %d-by-%d array with complex elements.', Axis.count, Nh, Nv);
			
			this.Jh = H_cell{v};
			this.Jv = -H_cell{h};
			
			this.Ft2d = cell(GK.count, Axis.count);
			for w = Axis.elems
				this.Ft2d{GK.prim, w} = array2scalar(H_cell{w}, PhysQ.H, this.grid2d, w, GK.prim, osc, this.intercept);
			end
			for w = Axis.elems
				this.Ft2d{GK.dual, w} = array2scalar(E_cell{w}, PhysQ.E, this.grid2d, w, GK.dual, osc, this.intercept);
			end
		end
		
		function [index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
			assert(~isempty(this.Jh) && ~isempty(this.Jv), '"Jh" and "Jv" are not set in this ModalSrc.');
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
				
				KA_curr = abs(this.Jh) .* dVh + abs(this.Jv) .* dVv;
				KA_curr = sum(KA_curr(:));
				norm_factor = this.KA/KA_curr;  % normalization factor
				
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

