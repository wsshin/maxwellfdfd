%% EllipticCylinder
% Concrete subclass of <GenenicCylinder.html |GenericCylinder|> representing a
% cylinder with an elliptical cross section.

%%% Description
% |EllipticCylinder| represents the shape of an ellicptic cylinder.  The major
% and minor axes of the ellipse as well as the axis of the cylinder should be
% aligned with the axes of the Cartesian coordinate system.

%%% Construction
%  shape = EllipticCylinder(normal_axis, height, center, semiaxes)
%  shape = EllipticCylinder(normal_axis, height, center, semiaxes, dl_max)
% 
% *Input Arguments*
%
% * |normal_axis|: axis of the cylinder.  It should be one of |Axis.x|,
% |Axis.y|, |Axis.z|.
% * |height|: size of the cylinder along its axis.
% * |center|: center of the cylinder in the format of |[x y z]|.
% * |semiaxes|: semiaxes of the ellipse in the format of |[a b]|.  If
% |normal_axis == Axis.y|, |a| is the semiaxis in the z-axis and |b| is the
% semiaxis in the x-axis.
% * |dl_max|: maximum grid size allowed in the cylinder.  It can be either |[dx
% dy dz]| or a single real number |dl| for |dx = dy = dz|.  If unassigned,
% |dl_max = Inf| is used.

%%% Example
%   % Create an instance of EllipticCylinder.
%   shape = EllipticCylinder(Axis.z, 100, [0 0 50], [100 50]);
%
%   % Use the constructed shape in maxwell_run().
%   [E, H] = maxwell_run({INITIAL ARGUMENTS}, 'OBJ', {'vacuum', 'none', 1.0}, shape, {REMAINING ARGUMENTS});

%%% See Also
% <CircularCylinder.html |CircularCylinder|>, <CircularShellCylinder.html
% |CircularShellCylinder|>, <SectoralCylinder.html |SectoralCylinder|>,
% <PolyognalCylinder.html |PolygonalCylinder|>, <Shape.html |Shape|>,
% <maxwell_run.html |maxwell_run|>

classdef EllipticCylinder < GenericCylinder

	methods
        function this = EllipticCylinder(normal_axis, height, center, semiaxes, dl_max)
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(height, 'real') && height > 0, '"height" should be positive.');
			chkarg(istypesizeof(center, 'real', [1, Axis.count]), ...
				'"center" should be length-%d row vector with real elements.', Axis.count);
			chkarg(istypesizeof(semiaxes, 'real', [1, Dir.count]) && all(semiaxes > 0), ...
				'"semiaxes" should be length-%d row vector with positive elements.', Axis.count);

			[h, v, n] = cycle(normal_axis);  % h-axis, v-axis, n-axis
			s = NaN(1, Axis.count);  % semisides
			s(h) = semiaxes(Dir.h);
			s(v) = semiaxes(Dir.v);
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
				
				level = zeros(size(p));
				for d = Dir.elems
					level = level + ((loc{d}-c(d)) ./ semiaxes(d)).^2;
				end
				level = 1 - sqrt(level);
			end
			
			lprim = cell(1, Axis.count);
			for w = Axis.elems
				lprim{w} = bound(w,:);
			end
			
			if nargin < 5  % no dl_max
				super_args = {normal_axis, @lsf2d, lprim};
			else
				super_args = {normal_axis, @lsf2d, lprim, dl_max};
			end
			
			this = this@GenericCylinder(super_args{:});
		end
	end
end

