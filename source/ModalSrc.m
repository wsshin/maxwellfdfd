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
%  src = ModalSrc(normal_axis, intercept, opts)
%  src = ModalSrc(normal_axis, intercept, opts, KA)
% 
% *Input Arguments*
%
% * |normal_axis|: direction normal to the plane over which the electric dipoles
% distribute.  It is also the axis of the waveguide.  It should be one of
% |Axis.x|, |Axis.y|, |Axis.z|.
% * |intercept|: location of the plane in the |normal_axis| direction.
% * |opts|: options structure that controls the behavior of the mode solver.
% See below for more details.
% * |KA|: integral of the norm of the surface current density.  It indicates the
% strength of the mode amplitude.  If unassigned, the deault value |KA = 1| is
% used.
%
% *Options Structure* 
%
% |opts| should have a string property |opts.clue| that indicates the type of
% the clue that the mode solver utilizes to find an appropriate waveguide mode.
% The clue can be either the order number of the waveguide mode, or the
% effective refractive index of the waveguide mode.
%
% * |opts.clue = 'order'| commands the mode solver to find the |n|th-order mode,
% where |n| is the order number you should provide as |opts.order = n|.  For
% example, |opts| with |opts.clue = 'order'| and |opts.order = 1| commands the
% mode solver to find the fundamental mode.
% * |opts.clue = 'guess'| commands the mode solver to find the mode with the
% effective index close to a "guess".  The guess of the effective index should
% be provided as |opts.neff = neff_guess|.  For example, |opts| with |opts.clue
% = 'guess'| and |opts.neff = 2.0 - 0.01i| commands the mode solver to find the
% mode whose effective index is close to |2.0 - 0.01i|.
%
% In general, |opts.clue = 'guess'| finds the mode much faster than |opts.clue =
% 'order'|.

%%% Note
% In the finite-difference grid, |ModalSrc| excites dipoles at the _E_-field
% points.  This poses a condition on |intercept| argument in the constructor:
% |intercept| should be at a dual grid point in the |normal_axis| direction.
% Therefore, make sure that |intercept| does not overlap with the locations of
% the vertices of <Shape.html |Shape|> in the |normal_axis| direction; otherwise
% dynamic grid generation in <moxwell_run.html |maxwell_run|> will fail.

%%% Example
%   % Create an instance of PointSrc.
%   modeopts.clue = 'order';
%   modeopts.order = 1;
%   src =  ModalSrc(Axis.y, -1000, modeopts);  % y = -1000 should not be primary grid point
%
%   % Use the constructed src in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'SRCJ', src);

%%% See Also
% <Plane.html |PlaneSrc|>, <TFSFPlaneSrc.html |TFSFPlaneSrc|>, <ModalSrc.html
% |ModalSrc|>, <maxwell_run.html |maxwell_run|>

classdef ModalSrc < Source
	
	properties (SetAccess = immutable)
		normal_axis  % plane normal axis: one of Axis.x, Axis.y, Axis.z
		intercept  % intercept between plane and normal axis
		KA  % surface integral of surface current density K; for z-normal plane, volume integral of (|Jx|+|Jy|)
		opts  % options to control method to calculate mode
	end
	
	properties (SetAccess = private)
		grid2d  % instance of Grid2d
		osc  % instance of Oscillation
		E2d  % {Ex2d Ey2d Ez2d}: cell array of Scalar2d for E on this plane
		H2d  % {Hx2d Hy2d Hz2d}: cell array of Scalar2d for H on this plane
		JMh  % 2D array J or M in horizontal direction on this plane: e.g., Jx for normal == z
		JMv  % 2D array J or M in vertical direction on this plane: e.g., Jy for normal == z
		neff  % effective n
	end
	
	% Properties for dispersion relation and propagation length
	properties (Dependent, SetAccess = immutable)
		ispreped  % true if this ModalSrc is prepared; false otherwise
		beta_r  % real part of complex wavevector
		Lp  % propagation length
	end
	
	methods
		function this = ModalSrc(normal_axis, intercept, opts, KA)
			chkarg(istypesizeof(normal_axis, 'Axis'), ...
				'"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			
			if nargin < 3  % no opts
				opts.clue = 'order';
				opts.order = 1;  % fundamental mode
			end
			chkarg(istypesizeof(opts, 'struct'), '"opts" should be structure.');

			chkarg(isequal(opts.clue, 'guess') || isequal(opts.clue, 'order'), ...
				'"opts.clue" should be ''guess'' or ''order'' (string).');
			if isequal(opts.clue, 'guess')
				chkarg(isfield(opts, 'neff') && istypesizeof(opts.neff, 'complex'), ...
					'"opts.neff" should be complex.');
				% if opts.clue == 'guess', opts can have additional field
				% opts.H2d.
			else  % opts.clue == 'order'
				if ~isfield(opts, 'order')
					opts.order = 1;  % fundamental mode
				end
				chkarg(istypesizeof(opts.order, 'int') && opts.order > 0, ...
					'"opts.order" should be positive integer.');
			end
			
			if nargin < 4  % no KA
				KA = 1.0;
			end
			chkarg(istypesizeof(KA, 'real'), '"KA" should be real.');
			
			lgrid = cell(1, Axis.count);
			laltgrid = cell(1, Axis.count);
			lgrid{normal_axis} = intercept;
			plane = Plane(normal_axis, intercept);
			this = this@Source(lgrid, laltgrid, plane);
			
			this.normal_axis = normal_axis;
			this.intercept = intercept;
			this.opts = opts;
			this.KA = KA;
			this.JMh = [];
			this.JMv = [];
		end
		
		function beta_r = get.beta_r(this)
			beta = 2*pi*this.neff / this.osc.in_L0();
			beta_r = real(beta);
		end
		
		function Lp = get.Lp(this)
			beta = 2*pi*this.neff / this.osc.in_L0();
			Lp = -1/imag(beta);
		end
		
		function ispreped = get.ispreped(this)
			ispreped = ~isempty(this.JMh) && ~isempty(this.JMh);
		end
		
		function setEH(this, neff, osc, E_cell, H_cell, ge, grid3d)
			chkarg(istypesizeof(neff, 'complex'), '"neff" should be complex.');
			this.neff = neff;
			
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			this.osc = osc;
			
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			this.grid2d = Grid2d(grid3d, this.normal_axis);

			chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');

			Nh = this.grid2d.N(Dir.h);
			Nv = this.grid2d.N(Dir.v);
			h = this.grid2d.axis(Dir.h);
			v = this.grid2d.axis(Dir.v);
			
			assert(istypesizeof(E_cell, 'complexcell', [1 Axis.count], [Nh Nv]), ...
				'"E_cell" should be length-%d row cell array whose each element is %d-by-%d array with complex elements.', Axis.count, Nh, Nv);
			assert(istypesizeof(H_cell, 'complexcell', [1 Axis.count], [Nh Nv]), ...
				'"H_cell" should be length-%d row cell array whole each element is %d-by-%d array with complex elements.', Axis.count, Nh, Nv);
			
			if this.gt == ge  % source is J
				this.JMh = H_cell{v};
				this.JMv = -H_cell{h};
			else  % source is M
				this.JMh = E_cell{v};
				this.JMv = -E_cell{h};
			end
			
			this.E2d = cell(1, Axis.count);
			this.H2d = cell(1, Axis.count);
			for w = Axis.elems
				gt = ge;  % grid type for E-field
				gt_array = gt(ones(1, Axis.count));
				gt_array(w) = alter(gt);
				this.E2d{w} = array2scalar(E_cell{w}, PhysQ.E, this.grid2d, w, FT.e, gt_array(this.grid2d.axis), osc, this.intercept);

				gt = alter(ge);  % grid type for H-field
				gt_array = gt(ones(1, Axis.count));
				gt_array(w) = alter(gt);
				this.H2d{w} = array2scalar(H_cell{w}, PhysQ.H, this.grid2d, w, FT.h, gt_array(this.grid2d.axis), osc, this.intercept);
			end
		end
		
		function [index_cell, JMw_patch] = generate_kernel(this, w_axis, grid3d)
			assert(~isempty(this.JMh) && ~isempty(this.JMv), '"JMh" and "JMv" are not set in this ModalSrc.');
			if w_axis == this.normal_axis
				JMw_patch = [];
				index_cell = cell(1, Axis.count);
			else
				g2d = Grid2d(grid3d, this.normal_axis);
				assert(isequal(g2d, this.grid2d), ...
					'%s-normal cross section of "grid3d" is different from the one set with JMh and JMv.', char(this.normal_axis));

				h = this.grid2d.axis(Dir.h);
				v = this.grid2d.axis(Dir.v);
				n = this.normal_axis;
				
				g = this.gt;
				ind_n = ind_for_loc(this.intercept, n, g, grid3d);
				
				% Set index_cell.
				index_cell = {':', ':', ':'};
				index_cell{n} = ind_n;
				
				% Set Jw_patch.
				dn = grid3d.dl{n,g}(ind_n);
				dha = grid3d.dl{h, alter(this.gt)};
				dhg = grid3d.dl{h, this.gt};
				dva = grid3d.dl{v, alter(this.gt)};
				dvg = grid3d.dl{v, this.gt};
				
				dVh = dn .* (dha.' * dvg);
				dVv = dn .* (dhg.' * dva);
				
				KA_curr = abs(this.JMh) .* dVh + abs(this.JMv) .* dVv;
				KA_curr = sum(KA_curr(:));  % KA_curr is real
				JMhv = [this.JMh(:); this.JMv(:)];
				[~, i_pf] = max(abs(JMhv));
				phasefactor = JMhv(i_pf)/abs(JMhv(i_pf));
				KA_curr = KA_curr * phasefactor;  % KA_curr is complex
				norm_factor = this.KA/KA_curr;  % normalization factor
				
				if w_axis == h
					JMw_patch = norm_factor .* this.JMh;
				else
					assert(w_axis == v);
					JMw_patch = norm_factor .* this.JMv;
				end
				
% 				if h > v  % h == Axis.z and v == Axis.x if this.normal_axis == Axis.y
% 					Jw_patch = permute(Jw_patch, int([Dir.v, Dir.h]));
% 				end
				JMw_patch = ipermute(JMw_patch, int([h v n]));
			end
		end
	end
	
end

