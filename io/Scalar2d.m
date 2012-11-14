classdef Scalar2d
    % SCALAR2D represents a 2D arary of scalar values.  It can represents a single
    % directional component of a 3D vector field (e.g. Ex of E).
    
    properties (SetAccess = immutable)
		% numerical properties
        array  % 2D array
        grid2d  % instance fo Grid2d
        osc  % instance of Oscillation
		intercept  % intercept of the plane of this 2D data in the normal axis

		% extra properties
		physQcell  % {PhysQ, int; PhysQ, int; ...}: representation of physical quanity
		unitvalue  % unit value of physical quantity
		name  % name of this quantity: "E_x", "H_y", etc. in LaTeX equation format
	end
	            
    methods
        function this = Scalar2d(array, grid2d, osc, physQcell, name, intercept)
			chkarg(istypesizeof(grid2d, 'Grid2d'), '"grid2d" should be instance of Grid2d.');
			this.grid2d = grid2d;
			
			N = grid2d.N + 1;
			chkarg(istypesizeof(array, 'complex', N), ...
				'"array" should be %d-by-%d array with complex elements.', N(Dir.h), N(Dir.v));
			this.array = array;
			
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			this.osc = osc;
			
			if nargin < 4  % no physQ
				physQcell = {PhysQ.arbitrary, 1};
			end
			this.physQcell = physQcell;
			this.unitvalue = osc.unit.value(this.physQcell);
			
			if nargin < 5  % no name
				name = '';
			end
            this.name = name;
			
			if nargin < 6  % no intercept
				intercept = NaN;
			end
			this.intercept = intercept;
		end
		
		function l_cell = l_data(this, withpml)
			% Return the locations where data are evaluated.
			l_cell = this.grid2d.lplot(GK.prim, withpml);
		end
		
		function l_cell = l_pixelbound(this, withpml)
			% Returen the locations of boundaries of pixels.
			l_cell = this.grid2d.lplot(GK.dual, withpml);
		end
		
		function [C, X, Y] = data_for_pcolor(this, withinterp, withpml)
			chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
			chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');
			
			C = this.array;
			if ~withpml
				Npml = this.grid2d.Npml;
				C = C(1+Npml(Dir.h,Sign.n):end-Npml(Dir.h,Sign.p), ...
					1+Npml(Dir.v,Sign.n):end-Npml(Dir.v,Sign.p));
			end
			
			if ~withinterp
				% Pad the end boundaries.  The padded values are not drawn, so
				% they don't have to be accurate.
				C(end+1, :) = C(end, :);
				C(:, end+1) = C(:, end);
			end

			% Permute dimensions to be compatible with pcolor().
			C = permute(C, int([Dir.v, Dir.h]));
% 			C = C .* this.unitvalue;
			
			if nargout >= 2  % X, Y
				if withinterp
					l = this.l_data(withpml);
				else
					l = this.l_pixelbound(withpml);
				end
				[X, Y] = meshgrid(l{:});
			end
		end
			
		function val = value(this, point)
			chkarg(istypesizeof(point, 'real', [0, Dir.count]), ...
				'"point" should be 2D array with %d columns.', Dir.count);
			chkarg(this.grid2d.contains(point), '"point" should be inside grid.');
			[C, X, Y] = this.data_for_pcolor(true, true);
			val = interp2(X, Y, C, point(Dir.h), point(Dir.v));
		end
	end		
end
