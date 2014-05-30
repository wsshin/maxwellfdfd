%% PlaneSrc
% Concrete subclass of <Source.html |Source|> representing a constant electric
% dipole distribution over an entire plane.

%%% Description
% |PlaneSrc| is to generate a plane wave in a homogeneous medium.  It supports
% oblique incidence, i.e., emission in the direction different from the normal
% direction of the plane.

%%% Construction
%  src = PlaneSrc(normal_axis, intercept, polarization)
%  src = PlaneSrc(normal_axis, intercept, polarization, K)
%  src = PlaneSrc(normal_axis, intercept, polarization, K, k_Bloch)
%  src = PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)
% 
% *Input Arguments*
%
% * |normal_axis|: direction normal to the plane.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |intercept|: location of the plane in the |normal_axis| direction.
% * |polarization|: direction of the dipoles distributed on the plane.  It
% should be either one of |Axis.x|, |Axis.y|, |Axis.z|, or an angle in radian.
% When it is an instance of <Axis.html |Axis|>, it should be orthogonal to
% |normal_axis|.  When it is an angle in radian, it is measured in the 2D
% Cartesian coordinate system normal to |normal_axis|.  For example, if
% |normal_axis == Axis.y|, then the angle is measured in the zx plane, and
% therefore it is measured in the counterclockwise direction from the z-axis in
% the zx plane.
% * |K|: amplitude of the surface current density that the dipoles drive.
% * |kt_Bloch|: in-plane Bloch wavevector.  If unassigned, |kt_Bloch = [0 0]| is
% used.
% * |theta|: oblique incidence angle measured from |normal_axis|.  |abs(theta)|
% should not exceed |pi/2| because the waves are emitted from the both sides of
% the plane.
% * |wvlen|: wavelength of the plane wave in the background medium.  It is
% required to set up the Bloch boundary condition for oblique incidence.
% |wvlen| is not the vacuum wavelength used in the frequency-domain Maxwell's
% equations, but the wavelength in the medium where this |PlaneSrc| is located.

%%% Note
% In the finite-difference grid, |PlaneSrc| excites dipoles at the _E_-field
% points.  This poses a condition on |intercept| argument in the constructor:
% |intercept| should be at a dual grid point in the |normal_axis| direction.
% Therefore, make sure that |intercept| does not overlap with the locations of
% the vertices of <Shape.html |Shape|> in the |normal_axis| direction; otherwise
% dynamic grid generation in <moxwell_run.html |maxwell_run|> will fail.

%%% Example
%   % Create an instance of PointSrc.
%   src =  PlaneSrc(Axis.y, 0, Axis.z);  % y = 0 should not be primary grid point
%
%   % Use the constructed src in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'SRCJ', src);

%%% See Also
% <PointSrc.html |PointSrc|>, <PointSrcM.html |PointSrcM|>, <TFSFPlaneSrc.html
% |TFSFPlaneSrc|>, <ModalSrc.html |ModalSrc|>, <maxwell_run.html |maxwell_run|>

classdef PlaneSrc < Source & WithBloch
	
	properties (SetAccess = immutable)
		normal_axis  % plane normal axis: one of Axis.x, Axis.y, Axis.z
		intercept  % intercept between plane and normal axis
		phi  % angle of polarization with respect to first Cartesian direction in this plane
		kBloch  % Bloch k-vector [kx, ky, kz]; for PhC mode solver, need to change the code to assign out-of-plane k component (but mode solver does not need source, so should be set directly to Grid)
		K  % surface current density
	end
	
	methods
		function this = PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)  % refractive index would be better rather than wvlen
			chkarg(istypesizeof(normal_axis, 'Axis'), ...
				'"normal_axis" should be instance of Axis.');
			[p, q] = cycle(normal_axis);

			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');			
			chkarg(istypesizeof(polarization, 'Axis') || istypesizeof(polarization, 'real'), ...
				'"polarization" should be instance of Axis or angle in radian.');
			phi_ = polarization;
			if istypesizeof(polarization, 'Axis')
				chkarg(polarization ~= normal_axis, '"polarization" should be orthogonal to "normal_axis".');
				if polarization == p
					phi_ = 0;
				else
					assert(polarization == q);
					phi_ = pi/2;
				end
			end
			
			if nargin < 4  % no K
				K = 1.0;
			end
			chkarg(istypesizeof(K, 'complex'), '"K" should be complex.');

			if nargin < 5  % no theta (or k_Bloch)
				kt_Bloch = [0 0];
			elseif nargin < 6  % no wvlen
				kt_Bloch = theta;  % theta is in fact kt_Bloch
				chkarg(istypesizeof(kt_Bloch, 'real', [1 Dir.count]), ...
					'"kt_Bloch" should be length-%d row vector with real elements.', Dir.count);
			else  % with wvlen
				assert(nargin == 6);
				chkarg(istypesizeof(theta, 'real') && theta >= -pi/2 && theta <= pi/2, ...
					'"theta" should be polar angle in radian between -pi/2 and pi/2.');
				chkarg(istypesizeof(wvlen, 'real') && wvlen > 0, '"wvlen" should be positive.');
				kt_Bloch = NaN(1, Dir.count);
				kt = (2*pi/wvlen) * sin(theta);  % k tangential to plane
				kp = kt * cos(phi_ + pi/2);
				kq = kt * sin(phi_ + pi/2);
				kt_Bloch(Dir.h) = kp;
				kt_Bloch(Dir.v) = kq;
			end
			
			lgrid = cell(1, Axis.count);
			laltgrid = cell(1, Axis.count);
			lgrid{normal_axis} = intercept;
			plane = Plane(normal_axis, intercept);
			this = this@Source(lgrid, laltgrid, plane);
			
			this.normal_axis = normal_axis;
			this.intercept = intercept;
			this.phi = phi_;
			this.K = K;
			this.kBloch = [0 0 0];
			this.kBloch([p q]) = kt_Bloch;
		end
		
		function [index_cell, JMw_patch] = generate_kernel(this, w_axis, grid3d)
			grid3d.set_kBloch(this);
			if w_axis == this.normal_axis
				JMw_patch = [];
				index_cell = cell(1, Axis.count);
			else
				n = this.normal_axis;
				w = w_axis;
				v = setdiff(Axis.elems, [n, w]);
				
				g = this.gt;
				ind_n = ind_for_loc(this.intercept, n, g, grid3d);
				dn = grid3d.dl{n,g}(ind_n);
				J = this.K / dn;  % for t normal to both n and K, K*dt = (current through dn*dt) = J * (dn*dt)
				
				if w == cycle(this.normal_axis)
					Jw = J * cos(this.phi);
				else
					Jw = J * sin(this.phi);
				end
				
				if Jw == 0
					JMw_patch = [];
					index_cell = cell(1, Axis.count);
				else
					%  Set Jw_patch.
					lw = grid3d.l{w, alter(g)};
					lv = grid3d.l{v, g};
					kw = this.kBloch(w);
					kv = this.kBloch(v);

% 					if w < v
% 						Jw_patch = Jw .* (exp(-1i * (kw .* lw)).' * exp(-1i * (kv .* lv)));
% 					else
% 						assert(w > v);
% 						Jw_patch = Jw .* (exp(-1i * (kv .* lv)).' * exp(-1i * (kw .* lw)));
% 					end

					JMw_patch = Jw .* (exp(-1i * (kw .* lw)).' * exp(-1i * (kv .* lv)));
					JMw_patch = ipermute(JMw_patch, int([w v n]));

					% Set index_cell.
					index_cell = {':', ':', ':'};
					index_cell{n} = ind_n;
				end
			end
		end		
	end
end

