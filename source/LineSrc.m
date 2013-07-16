%% LineSrc
% Concrete subclass of <Source.html |Source|> representing a constant electric
% dipole distribution along a line.

%%% Description
% |LineSrc| is to generate a cylindrical wave in a homogeneous medium.  It can
% creates a cylindrical wave propagating not only in the radial direction, but
% also in the direction slanted toward the line axis.  It also supports the
% polarization normal to the line axis.

%%% Construction
%  src = LineSrc(axis, intercept)
%  src = LineSrc(axis, intercept, polarization)
%  src = LineSrc(axis, intercept, polarization, IorK)
%  src = LineSrc(axis, intercept, polarization, IorK, ka_Bloch)
%  src = LineSrc(axis, intercept, polarization, IorK, theta, wvlen)
% 
% *Input Arguments*
%
% * |axis|: direction of the line.  It should be one of |Axis.x|, |Axis.y|,
% |Axis.z|.
% * |intercept|: location of the line in the plane normal to the line.  For
% |axis == Axis.y|, it is in the format of |[z x]|.
% * |polarization|: direction of the dipoles distributed on the plane.  It
% should be one of |Axis.x|, |Axis.y|, |Axis.z|.  If unassigned, |polarization =
% axis| is used.
% * |IorK|: for |polarization == axis|, complex amplitude of the current that
% the dipoles drive; for |polarization ~= axis|, complex amplitude of the sheet
% current density that the dipoles drive.  If unassigned, |IorK = 1.0| is used.
% * |ka_Bloch|: Bloch wavevector in the axis of this line.  If unassigned,
% |ka_Bloch = 0| is used.
% * |theta|: oblique incidence angle.  |abs(theta)| should not exceed |pi/2|. If
% unassigned, |theta = 0.0| is used.
% * |wvlen|: wavelength of the cylindrical wave in the background medium.  It is
% required to set up the Bloch boundary condition for oblique incidence. |wvlen|
% is not the vacuum wavelength used in the frequency-domain Maxwell's equations,
% but the wavelength in the medium where this |LineSrc| is located.

%%% Note
% In the finite-difference grid, the dipoles in |LineSrc| are located at the
% _E_-field points.  This poses a condition on |intercept| argument in the
% constructor: if |polarization == axis|, both elements of |intercept| should be
% at dual grid points; if |polarization ~= axis|, the element of |intercept|
% corresponding to the |polarization| component should be at a primary grid
% point, and the other element of |intercept| should be at a dual grid point.
% Therefore, make sure that the elements of |intercept| at primary grid points
% do not overlap with the locations of the vertices of <Shape.html |Shape|>;
% otherwise dynamic grid generation in <moxwell_run.html |maxwell_run|> will
% fail.

%%% Example
%   % Create an instance of LineSrc.
%   src = LineSrc(Axis.z, [0 0]);
%
%   % Use the constructed src in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'SRCJ', src);

%%% See Also
% <PointSrc.html |PointSrc|>, <PlaneSrc.html |PlaneSrc|>, <maxwell_run.html
% |maxwell_run|>

classdef LineSrc < Source & WithBloch
	
	properties (SetAccess = immutable)
		axis  % axis of line: one of Axis.x, Axis.y, Axis.z
		intercept  % [h, v]: location of intercept in transverse plane
		polarization  % one of Axis.x, Axis.y, Axis.z
		kBloch  % Bloch k-vector [kx, ky, kz]
		IorK  % current I if polarization == axis; sheet current density K otherwise
	end
	
	methods
		function this = LineSrc(axis, intercept, polarization, IorK, theta, wvlen)
			chkarg(istypesizeof(axis, 'Axis'), '"axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real', [1 Dir.count]), ...
				'"intercept" should be length-%d row vector with real elements.', Dir.count);
			
			if nargin < 3  % no polarization
				polarization = axis;
			end
			chkarg(istypesizeof(polarization, 'Axis'), '"polarization" should be instance of Axis.');

			if nargin < 4  % no IorK
				IorK = 1.0;
			end
			chkarg(istypesizeof(IorK, 'complex'), '"IorK" should be complex.');

			if nargin < 5  % no theta (or ka_Bloch)
				ka_Bloch = 0;
			elseif nargin < 6  % no wvlen
				ka_Bloch = theta;  % theta is in fact ka_Bloch
				chkarg(istypesizeof(ka_Bloch, 'real'), '"ka_Bloch" should be real.');
			else  % with wvlen
				assert(nargin == 6);
				chkarg(istypesizeof(theta, 'real') && theta >= -pi/2 && theta <= pi/2, ...
					'"theta" should be polar angle in radian between -pi/2 and pi/2.');
				chkarg(istypesizeof(wvlen, 'real') && wvlen > 0, '"wvlen" should be positive.');
				ka_Bloch = (2*pi/wvlen) * sin(theta);
			end
		
			lgrid = cell(1, Axis.count);
			laltgrid = cell(1, Axis.count);
			[h, v] = cycle(axis);  % h, v ~= axis
			if polarization == h
				laltgrid{h} = intercept(Dir.h);
				lgrid{v} = intercept(Dir.v);
			elseif polarization == v
				lgrid{h} = intercept(Dir.h);
				laltgrid{v} = intercept(Dir.v);
			else  % polarization == axis
				lgrid{h} = intercept(Dir.h);
				lgrid{v} = intercept(Dir.v);
			end
				
			line = Line(axis, intercept);
			this = this@Source(lgrid, laltgrid, line);
			
			this.axis = axis;
			this.intercept = intercept;
			this.polarization = polarization;
			this.IorK = IorK;

			this.kBloch = [0 0 0];
			this.kBloch(this.axis) = ka_Bloch;
		end
				
		function [index_cell, JMw_patch] = generate_kernel(this, w_axis, grid3d)
			grid3d.set_kBloch(this);
			[h, v, a] = cycle(this.axis);  % h, v: axes in plane normal to this line
			p = this.polarization;
			if w_axis ~= p
				JMw_patch = [];
				index_cell = cell(1, Axis.count);
			else  % w_axis == p
				index_cell = {':', ':', ':'};
				if p == a
					% Set index_cell{h} and index_cell{v}.
					g = this.gt;
					axes = [h v];
					for d = Dir.elems
						ind = ind_for_loc(this.intercept(d), axes(d), g, grid3d);
						index_cell{axes(d)} = ind;
					end
					
					% Set Jw.
					dh = grid3d.dl{h,g}(index_cell{h});
					dv = grid3d.dl{v,g}(index_cell{v});
					I = this.IorK;
					JMw = I / (dh * dv);  % I = J * (area)
					ga = alter(this.gt);
				else  % p ~= a
					% Set index_cell{h}.
					if h == p
						g = alter(this.gt);
					else
						g = this.gt;
					end
					ind = ind_for_loc(this.intercept(Dir.h), h, g, grid3d);
					index_cell{h} = ind;
					
					% Set index_cell{v}
					g = alter(g);
					ind = ind_for_loc(this.intercept(Dir.v), v, g, grid3d);
					index_cell{v} = ind;
					
					% Set Jw.
					q = setdiff([h v], p);
					dq = grid3d.dl{q, this.gt}(index_cell{q});
					K = this.IorK;
					JMw = K / dq;
					ga = this.gt;
				end
				
				% Set Jw_patch.
				la = grid3d.l{a, ga};
				ka = this.kBloch(a);

				JMw_patch = JMw .* exp(-1i * (ka .* la)).';
				JMw_patch = ipermute(JMw_patch, int([a h v]));
			end
		end		
	end
end
