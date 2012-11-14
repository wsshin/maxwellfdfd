classdef Scalar3d
    % SCALAR3D represents a 3D arary of scalar values.  It can represents a single
    % directional component of a 3D vector field (e.g. Ex of E).
	% The fields are always evaluated at the primary grid points (grid kind =
	% GK.prim)
    
    properties (SetAccess = immutable)
		% numerical properties
        array  % 3D array
        grid3d  % instance fo Grid3d
        osc  % instance of Oscillation

		% extra properties
		physQcell  % {PhysQ, int; PhysQ, int; ...}: representation of physical quanity
		unitvalue  % unit value of physical quantity
		name  % name of this quantity: "E_x", "H_y", etc. in LaTeX equation format
	end
	            
    methods
        function this = Scalar3d(array, grid3d, osc, physQcell, name)
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			this.grid3d = grid3d;
			
			N = grid3d.N + 1;
			chkarg(istypesizeof(array, 'complex', N), ...
				'"array" should be %d-by-%d-by-%d array with complex elements.', N(Axis.x), N(Axis.y), N(Axis.z));
			this.array = array;
			
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			this.osc = osc;
			
			if nargin < 4  % no physQ
				this.physQcell = {PhysQ.arbitrary, 1};
			end
			this.physQcell = physQcell;
			this.unitvalue = osc.unit.value(this.physQcell);
			
			if nargin < 5  % no name
				name = '';
			end
            this.name = name;
		end
		
		function l_cell = l_data(this, withpml)
			% Return the locations where data are evaluated.
			l_cell = this.grid3d.lplot(GK.prim, withpml);
		end
		
		function l_cell = l_voxelbound(this, withpml)
			% Returen the locations of boundaries of pixels.
			l_cell = this.grid3d.lplot(GK.dual, withpml);
		end
		
		function [V, X, Y, Z] = data_for_slice(this, withinterp, withpml)
			chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
			chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');
			
			V = this.array;
			if ~withpml
				Npml = this.grid3d.Npml;
				V = V(1+Npml(Axis.x,Sign.n):end-Npml(Axis.x,Sign.p), ...
					1+Npml(Axis.y,Sign.n):end-Npml(Axis.y,Sign.p), ...
					1+Npml(Axis.z,Sign.n):end-Npml(Axis.z,Sign.p));
			end
			
			if ~withinterp
				% Pad the end boundaries.  The padded values are not drawn, so
				% they don't have to be accurate.
				V(end+1, :, :) = V(end, :, :);
				V(:, end+1, :) = V(:, end, :);
				V(:, :, end+1) = V(:, :, end);
			end

			% Permute dimensions to be compatible with meshgrid().
			V = permute(V, int([Axis.y, Axis.x, Axis.z]));
% 			V = V .* this.unitvalue;
			
			if nargout >= 2  % X, Y, Z
				if withinterp
					l = this.l_data(withpml);
				else
					l = this.l_voxelbound(withpml);
				end
				[X, Y, Z] = meshgrid(l{:});
			end
		end
			
		function val = value(this, point)
			chkarg(istypesizeof(point, 'real', [0, Axis.count]), ...
				'"point" should be 2D array with %d columns.', Axis.count);
			chkarg(this.grid3d.contains(point), '"point" should be inside grid.');
			[V, X, Y, Z] = this.data_for_slice(true, true);
			val = interp3(X, Y, Z, V, point(Axis.x), point(Axis.y), point(Axis.z));
		end
	end		
end
