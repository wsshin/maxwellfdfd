%% PointSrc
% Concrete subclass of <Source.html |Source|> representing an electric point
% dipole source.

%%% Description
% |PointSrc| is the most basic source type.  It places an oscillating electric
% point dipole at the location given in the constructor.

%%% Construction
%  src = PointSrc(polarization_axis, location)
%  src = PointSrc(polarization_axis, location, I)
% 
% *Input Arguments*
%
% * |polarization_axis|: direction of the dipole.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |location|: location of the dipole in the form of |[x, y, z]|.
% * |I|: amplitude of the current that the dipole drives.

%%% Note
% In the finite-difference grid, |PointSrc| is located at one of the _E_-field
% points.  This poses a condition on |location| argument in the constructor: the
% location in the directions normal to the dipole polarization should be at dual
% grid points, whereas the location in the direction along the dipole
% polarization should be at a primary grid point.  Therefore, make sure that the
% location of the dipole does not overlap with the locations of the vertices of
% <Shape.html |Shape|> in the directions normal to the dipole polarization;
% otherwise dynamic grid generation in <moxwell_run.html |maxwell_run|> will
% fail.

%%% Example
%   % Create an instance of PointSrc.
%   src = PointSrc(Axis.z, [0 0 0]);  % z = 0 should not be dual grid point
%
%   % Use the constructed src in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'SRC', src);

%%% See Also
% <PointSrcM.html |PointSrcM|>, <PlaneSrc.html |PlaneSrc|>, <maxwell_run.html
% |maxwell_run|>

classdef PointSrc < Source
	
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
					iv = ind_for_loc(l, v, g, grid3d);
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
