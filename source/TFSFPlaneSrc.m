classdef TFSFPlaneSrc < Source
	% TFSFSrc is a source used for the TF/SF method.
		
	properties (SetAccess = immutable)
		khat  % [nx, ny, nz]: unit vector along k-vector of plane wave
		Ehat  % [px, py, pz]: unit vector along polarization direction
		E0  % complex amplitude of plane wave
	end
	
	properties (SetAccess = private)
		N_bg  % complex refractive index of background medium
		grid3d  % Grid3d instance used to set up E and J
		J  % {Jx_array, Jy_array, Jz_array}: J generating plane wave in TF/SF box
	end
	
	methods
		function this = TFSFPlaneSrc(bound, propagation_axis, polarization_axis, E0)
			chkarg(istypesizeof(bound, 'real', [Axis.count, Sign.count]), ...
				'"bound" should be %d-by-%d array of real elements.', Axis.count, Sign.count);
			chkarg(istypesizeof(propagation_axis, 'Axis') || istypesizeof(propagation_axis, 'real', [1, Axis.count]), ...
				'"propagation_axis" should be instance of Axis or length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(polarization_axis, 'Axis') || istypesizeof(polarization_axis, 'real', [1, Axis.count]), ...
				'"polarization_axis" should be instance of Axis or length-%d row vector with real elements.', Axis.count);
			
			if nargin < 4
				E0 = 1;
			end
			
			l = cell(Axis.count, GK.count);
			for w = Axis.elems
				l{w,GK.prim} = bound(w,:);
			end
			this = this@Source(l, Box(bound));

			this.khat = propagation_axis;
			if istypesizeof(propagation_axis, 'Axis')
				this.khat = zeros(1, Axis.count);
				this.khat(propagation_axis) = 1;
			end
			this.khat = this.khat / norm(this.khat);
			
			this.Ehat = polarization_axis;
			if istypesizeof(polarization_axis, 'Axis')
				this.Ehat = zeros(1, Axis.count);
				this.Ehat(polarization_axis) = 1;
			end
			this.Ehat = this.Ehat / norm(this.Ehat);
			
			chkarg(this.khat * this.Ehat.' < 1e-12, '"normal_axis" and "polarization_axis" should be orthogonal.');
			
			this.E0 = E0;
			this.N_bg = NaN;
			this.grid3d = Grid3d.empty();
		end
		
		function set_bg_material(this, material)
			chkarg(istypesizeof(material, 'Material'), '"mat_bg" should be instance of Material.');
			this.N_bg = sqrt(material.eps * material.mu);
		end
		
		function E_cell = create_incidentE(this, osc, grid3d)
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			this.grid3d = grid3d;

			E_cell = cell(1, Axis.count);
			l = cell(1, Axis.count);
			bound = this.shape.bound;
			ind = cell(1, Axis.count);
			k = this.khat * 2*pi * this.N_bg / osc.in_L0();  % k*n = (2*pi*n/lambda)
			for w = Axis.elems
				E_cell{w} = zeros(grid3d.N);
				
				[p, q] = cycle(w);
				l{p} = grid3d.l{p, GK.dual};
				l{q} = grid3d.l{q, GK.dual};
				l{w} = grid3d.l{w, GK.prim};
				
				for v = Axis.elems
					ind_n = find(l{v} >= bound(v,Sign.n), 1, 'first');
					ind_p = find(l{v} <= bound(v,Sign.p), 1, 'last');
					ind{v} = ind_n:ind_p;
				end
				
				[X, Y, Z] = ndgrid(l{Axis.x}(ind{Axis.x}), l{Axis.y}(ind{Axis.y}), l{Axis.z}(ind{Axis.z}));
				E_cell{w}(ind{:}) = (this.E0 * this.Ehat(w)) * exp(1i * (k(Axis.x)*X + k(Axis.y)*Y + k(Axis.z)*Z));
			end
		end
		
		function setJ(this, J_cell, grid3d)
			chkarg(isequal(grid3d, this.grid3d), 'same "grid3d" should be used as in create_incidentE().');
			chkarg(istypesizeof(J_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
				'"J_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
				Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
			
			this.J = J_cell;
		end
		
		function [index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
			assert(~isnan(this.N_bg), '"mat_bg" is not set in this TFSFPlaneSrc.');
			
			index_cell = cell(1, Axis.count);
			for w = Axis.elems
				index_cell{w} = 1:grid3d.N(w);
			end
			Jw_patch = this.J{w_axis};
		end
	end
end
