%% CircularShellCylinder
% Concrete subclass of <GenenicCylinder.html |GenericCylinder|> representing a
% cylinder with a circular shell as the cross section (the cylinder is a hollow
% pipe).

%%% Description
% |CircularShellCylinder| represents the shape of a hollow cylindrical pipe with
% circular cross sections.  The axis of the cylinder should be aligned with an
% axis of the Cartesian coordinate system.

%%% Construction
%  shape = CircularShellCylinder(normal_axis, height, center, r1, r2)
%  shape = CircularShellCylinder(normal_axis, height, center, r1, r2, dl_max)
% 
% *Input Arguments*
%
% * |normal_axis|: axis of the cylinder.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |height|: size of the cylinder along its axis.
% * |center|: center of the cylinder in the format of |[x y z]|.
% * |r1|, |r2|: radii of the inner and outer circles (or vice versa).
% * |dl_max|: maximum grid size allowed in the cylinder.  It can be either |[dx
% dy dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned,
% |dl_max = Inf| is used.

%%% Example
%   % Create an instance of CircularShellCylinder.
%   shape = CircularShellCylinder(Axis.z, 100, [0 0 50], [100 50], 25, 50);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <CircularCylinder.html |CircularCylinder|>, <EllipticCylinder.html
% |EllipticCylinder|>, <SectoralCylinder.html |SectoralCylinder|>,
% <PolyognalCylinder.html |PolygonalCylinder|>, <Shape.html |Shape|>,
% <maxwell_run.html |maxwell_run|>

classdef CircularShellCylinder < GenericCylinder

	methods
        function this = CircularShellCylinder(normal_axis, height, center, r1, r2, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');
			chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
				'"center" should be length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(r1, 'real') && r1 > 0, '"r1" should be positive.');
			chkarg(istypesizeof(r2, 'real') && r2 > 0, '"r2" should be positive.');

			[h, v, n] = cycle(normal_axis);  % h-axis, v-axis, n-axis
			rs = sort([r1 r2]);
			r = rs(1);
			R = rs(2);
			
			s = NaN(1, Axis.count);  % semisides
			s(h) = R;
			s(v) = R;
			s(n) = height / 2;
			bound = [center - s; center + s];
			bound = bound.';
			
			% For c = center([h, v]), the level set function is
			% basically 1 - norm((rho-c) ./ semiaxes) but it is vectorized,
			% i.e., modified to handle rho = [p q] with column vectors p and q.
			function level = lsf2d(p, q)
				chkarg(istypeof(p, 'real'), '"p" should be array with real elements.');
				chkarg(istypeof(q, 'real'), '"q" should be array with real elements.');
				chkarg(isequal(size(p), size(q)), '"p" and "q" should have same size.');

				c = center([h, v]);				
				loc = {p, q};
				
				level1 = zeros(size(p));  % level set function for inside of R
				level2 = zeros(size(p));  % level set function for outside of r
				for d = Dir.elems
					level1 = level1 + ((loc{d}-c(d)) ./ R).^2;
					level2 = level2 + ((loc{d}-c(d)) ./ r).^2;
				end
				level1 = 1 - sqrt(level1);
				level2 = sqrt(level2) - 1;
				
				level = min(level1, level2);  % intersection of two regions
			end
			
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				lprim{w} = bound(w,:);
			end
			
			if nargin < 6  % no dl_max
				super_args = {normal_axis, @lsf2d, lprim};
			else
				super_args = {normal_axis, @lsf2d, lprim, dl_max};
			end
			
			this = this@GenericCylinder(super_args{:});
		end
	end
end

