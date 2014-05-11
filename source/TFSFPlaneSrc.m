%% TFSFPlaneSrc
% Concrete subclass of <Source.html |Source|> representing plane wave
% propagating in a total field box of the total-field/scattered-field (TF/SF)
% method.

%%% Description
% |TFSFPlaneSrc| is a source type for the TF/SF method.  It generates a plane
% wave inside a total field box.  The box should be empty except for scatterers.
% It supports an arbitrary propagation direction and polarization direction.

%%% Construction
%  src = TFSFPlaneSrc(bound, propagation_axis, polarization_axis)
%  src = TFSFPlaneSrc(bound, propagation_axis, polarization_axis, F0)
% 
% *Input Arguments*
%
% * |bound|: description of the box in the format of |[xmin xmax; ymin ymax;
% zmin zmax]|
% * |propagation_axis|: direction of the wavevector of the plane wave.  It
% can be either one of |Axis.x|, |Axis.y|, |Axis.z|, or |[kx ky kz]|.  The
% latter does not have to be normalized, i.e., |norm([kx ky kz])| does not have
% to be |1|.
% * |polarization_axis|: axis of the polarization of the _F_-field (the _E_- or
% _H_-field depending on whether the source is used as 'SRCJ' or 'SRCM') of the
% plane wave.  It can be either one of |Axis.x|, |Axis.y|, |Axis.z|, or |[Fx Fy
% Fz]|. The latter does not have to be normalized, i.e., |norm([Fx Fy Fz])| does
% not have to be |1|.
% * |F0|: complex amplitude of the plane wave measure at the origin. If not
% assigned, the default value |F0 = 1| is used.


%%% Methods
% * |set_bg_material(this, material)|: sets the background material that
% determines the wavelength of the plane wave.  When |TFSFPlaneSrc| is passed as
% a SRC paarmeter in <maxwell_run.html |maxwell_run|>, the material used as the
% default material of the simulation domain, i.e., the material assigned in the
% DOM parameter group, is taken automatically.
% * |F_cell = create_incidentF(this, osc, grid3d)|: returns |{Fx_array,
% Fy_array, Fz_array}| that has the _F_-field of the plane wave inside the total
% field box region.  When |TFSFPlaneSrc| is used as a SRC parameter in
% <maxwell_run.html |maxwell_run|>, <maxwell_run.html |maxwell_run|> calls this
% function to construct |J| or |M| that creates the _F_-field inside the total
% field box.

%%% Note
% |TFSFPlaneSrc| is similar to <PlaneSrc.html |PlaneSrc|>, but has a few
% differences. |PlaneSrc| requires a periodic boundary condition along the
% transverse direction, so it can be used only for a periodic structure.
% Therefore, to simulate a plane wave incident on a non-periodic, isolated
% scatterer, |TFSFPlaneSrc| sohould be used.  |PlaneSrc| supports an arbitrary
% propagation direction by utilizing the Bloch boundary condition internally,
% but it does not support an arbitrary polarization direction: the polarization
% direction is restricted to |Axis.x|, |Axis.y|, |Axis.z|.  On the other hand,
% |TFSFPlaneSrc| supports a fully arbitrary polarization direction.
%
% When passed as a SRC parameter of <maxwell_run.html |maxwell_run|>, the
% scatterer objects are described in the SOBJ parameter group and the background
% objects are described in the OBJ parameter group.  No background objects
% should be inside the total field box, because |TFSFPlaneSRC| assumes that the
% total field box is completely filled with the default material of the
% simulation domain that is described in the DOM parameter group of
% <maxwell_run.html |maxwell_run|>.

%%% Example
%   % Create an instance of TFSFPlaneSrc.
%   src =  TFSFPlaneSrc(Axis.y, 0, Axis.z);  % y = 0 should not be primary grid point
%
%   % Use the constructed src in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'SRCJ', src);

%%% See Also
% <PlaneSrc.html |PlaneSrc|>, <ModalSrc.html |ModalSrc|>, <maxwell_run.html |maxwell_run|>

classdef TFSFPlaneSrc < Source
		
	properties (SetAccess = immutable)
		khat  % [nx, ny, nz]: unit vector along k-vector of plane wave
		Fhat  % [px, py, pz]: unit vector along polarization direction
		F0  % complex amplitude of plane wave
	end
	
	properties (SetAccess = private)
		N_bg  % complex refractive index of background medium
		grid3d  % Grid3d instance used to set up F and J or M
		JM  % {JMx_array, JMy_array, JMz_array}: J or M generating plane wave in TF/SF box
	end
	
	methods
		function this = TFSFPlaneSrc(bound, propagation_axis, polarization_axis, F0)
			chkarg(istypesizeof(bound, 'real', [Axis.count, Sign.count]), ...
				'"bound" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
			chkarg(istypesizeof(propagation_axis, 'Axis') || istypesizeof(propagation_axis, 'real', [1, Axis.count]), ...
				'"propagation_axis" should be instance of Axis or length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(polarization_axis, 'Axis') || istypesizeof(polarization_axis, 'real', [1, Axis.count]), ...
				'"polarization_axis" should be instance of Axis or length-%d row vector with real elements.', Axis.count);
			
			if nargin < 4
				F0 = 1;
			end
			chkarg(istypesizeof(F0, 'complex'), '"F0" should be complex.');
			
			lgrid = cell(1, Axis.count);
			laltgrid = cell(1, Axis.count);
			for w = Axis.elems
				lgrid{w} = bound(w,:);
			end
			
			forceprim = true;
			this = this@Source(lgrid, laltgrid, Box(bound), forceprim);

			this.khat = propagation_axis;
			if istypesizeof(propagation_axis, 'Axis')
				this.khat = zeros(1, Axis.count);
				this.khat(propagation_axis) = 1;
			end
			this.khat = this.khat / norm(this.khat);
			
			this.Fhat = polarization_axis;
			if istypesizeof(polarization_axis, 'Axis')
				this.Fhat = zeros(1, Axis.count);
				this.Fhat(polarization_axis) = 1;
			end
			this.Fhat = this.Fhat / norm(this.Fhat);
			
			chkarg(this.khat * this.Fhat.' < 1e-12, '"normal_axis" and "polarization_axis" should be orthogonal.');
			
			this.F0 = F0;
			this.N_bg = NaN;
			this.grid3d = Grid3d.empty();
		end
		
		function set_bg_material(this, material)
			chkarg(istypesizeof(material, 'Material'), '"mat_bg" should be instance of Material.');
			this.N_bg = sqrt(material.eps * material.mu);
		end
		
		function F_cell = create_incidentF(this, osc, grid3d)
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			this.grid3d = grid3d;

			F_cell = cell(1, Axis.count);
			l = cell(1, Axis.count);
			bound = this.shape.bound;
			ind = cell(1, Axis.count);
			k = this.khat * 2*pi * this.N_bg / osc.in_L0();  % k = (2*pi*n/lambda)
			for w = Axis.elems
				F_cell{w} = zeros(grid3d.N);
				
				[p, q] = cycle(w);
				l{p} = grid3d.l{p, this.gt};
				l{q} = grid3d.l{q, this.gt};
				l{w} = grid3d.l{w, alter(this.gt)};
				
				for v = Axis.elems
					ind_n = find(l{v} >= bound(v,Sign.n), 1, 'first');
					ind_p = find(l{v} <= bound(v,Sign.p), 1, 'last');
					ind{v} = ind_n:ind_p;
				end
				
				[X, Y, Z] = ndgrid(l{Axis.x}(ind{Axis.x}), l{Axis.y}(ind{Axis.y}), l{Axis.z}(ind{Axis.z}));
				F_cell{w}(ind{:}) = (this.F0 * this.Fhat(w)) * exp(-1i * (k(Axis.x)*X + k(Axis.y)*Y + k(Axis.z)*Z));
			end
		end
		
		function setJM(this, JM_cell, grid3d)
			chkarg(isequal(grid3d, this.grid3d), 'same "grid3d" should be used as in create_incidentF().');
			chkarg(istypesizeof(JM_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
				'"JM_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements.', ...
				Axis.count, grid3d.N(Axis.x), grid3d.N(Axis.y), grid3d.N(Axis.z));
			
			this.JM = JM_cell;
		end
		
		function [index_cell, JMw_patch] = generate_kernel(this, w_axis, grid3d)
			assert(~isnan(this.N_bg), '"mat_bg" is not set in this TFSFPlaneSrc.');
			
			index_cell = {':', ':', ':'};
			JMw_patch = this.JM{w_axis};
		end
	end
end
