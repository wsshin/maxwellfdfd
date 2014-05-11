classdef Scalar3d
    % SCALAR3D represents a 3D arary of scalar values.  It can represents a single
    % directional component of a 3D vector field (e.g. Ex of E).
	% The fields are always evaluated at the primary grid points (grid kind =
	% GT.prim)
    
    properties (SetAccess = immutable)
		% numerical properties
        array  % 3D array
        grid3d  % instance of Grid3d
		gt_array  % [x_gridtype y_gridtype z_gridtype]: array of instances of GT
        osc  % instance of Oscillation

		% extra properties
		physQcell  % {PhysQ, int; PhysQ, int; ...}: representation of physical quanity
		unitvalue  % unit value of physical quantity
		name  % name of this quantity: "E_x", "H_y", etc. in LaTeX equation format
	end
	            
    methods
        function this = Scalar3d(array, grid3d, gt_array, osc, physQcell, name)
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			this.grid3d = grid3d;
			
			chkarg(istypesizeof(gt_array, 'GT', [1 Axis.count]), ...
				'"gt_array" should be length-%d row vector with GT as elements.', Axis.count);
			this.gt_array = gt_array;

			N = grid3d.N + int(gt_array);
			chkarg(istypesizeof(array, 'complex', N), ...
				'"array" should be %d-by-%d-by-%d array with complex elements.', N(Axis.x), N(Axis.y), N(Axis.z));
			this.array = array;			
			
			chkarg(istypesizeof(osc, 'Oscillation'), '"osc" should be instance of Oscillation.');
			this.osc = osc;
			
			if nargin < 5  % no physQ
				physQcell = PhysQ.arbitrary;
			end
			chkarg(istypesizeof(physQcell, 'PhysQ') || isphysQcell(physQcell), ...
				'"physQcell" should be instance of PhysQ, or cell array {PhysQ, int; PhysQ, int; ...}.');
			if istypesizeof(physQcell, 'PhysQ')
				physQcell = {physQcell, 1};
			end
			this.physQcell = physQcell;
			this.unitvalue = osc.unit.value(this.physQcell);
			
			if nargin < 6  % no name
				name = '';
			end
			chkarg(ischar(name), '"name" should be string.');
            this.name = name;
		end
		
		function l_cell = lplot(this, withinterp, withpml)
			% Return the locations where data are evaluated for plotting.  If
			% the data do not include the boundaries of the simulation domain (or
			% the PML interfaces for "withpml == false"), the boundary points
			% are added.
			l_cell = this.grid3d.lplot(this.gt_array, withinterp, withpml);
		end
		
		function l_cell = lvoxelbound(this, withpml)
			% Return the locations of boundaries of voxels drawn.  For data at
			% primary grid points, the voxel centers are in the simulation
			% domain including the boundary.  For data at dual grid points, the
			% voxel centers are in the simulation domain excluding the boundary.
			l_cell = this.grid3d.lvoxelbound(this.gt_array, withpml);
		end
		
		function [array, l_cell] = data_expanded(this)
			l_cell = this.grid3d.lall(Axis.elems + Axis.count*subsindex(this.gt_array));
			array = this.array;
		end
		
		function [array, l_cell] = data_original(this)
			l_cell = this.grid3d.l(Axis.elems + Axis.count*subsindex(this.gt_array));
			ind = cell(1, Axis.count);
			
			for w = Axis.elems
				if this.gt_array(w) == GT.prim
					ind{w} = 1:this.grid3d.N(w);
				else  % this.gt_array(w) == GT.dual
					ind{w} = 2:(this.grid3d.N(w)+1);
				end
			end
			array = this.array(ind{:});
		end
		
		function [V, X, Y, Z] = data_for_slice(this, withinterp, withpml)
			chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
			chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');
			
			lall = this.grid3d.lall(Axis.elems + Axis.count*subsindex(this.gt_array));
			[X, Y, Z] = ndgrid(lall{:});

			lplot = this.lplot(withinterp, withpml);
			[Xi, Yi, Zi] = ndgrid(lplot{:});
			V = interpn(X, Y, Z, this.array, Xi, Yi, Zi);
			V = permute(V, int([Axis.y, Axis.x, Axis.z]));  % to be compatible with slice()

			if ~withinterp
				% Pad the end boundaries.  The padded values are not drawn, so
				% they don't have to be accurate.
				V(end+1, :, :) = V(end, :, :);
				V(:, end+1, :) = V(:, end, :);
				V(:, :, end+1) = V(:, :, end);
			end
			
			if nargout >= 2  % X, Y, Z
				if withinterp
					l = this.lplot(withinterp, withpml);
				else
					l = this.lvoxelbound(withpml);
				end
				[X, Y, Z] = meshgrid(l{:});
			end
		end
		
		function val = value(this, point)
			chkarg(istypesizeof(point, 'real', [0, Axis.count]), ...
				'"point" should be 2D array with %d columns.', Axis.count);
			chkarg(this.grid3d.contains(point), '"point" should be inside grid.');
			
			lall = this.grid3d.lall(Axis.elems + Axis.count*subsindex(this.gt_array));
			ind = cell(1, Axis.count);
			for w = Axis.elems
				indw = ismembc2(point(w), lall{w});
				if indw == 0
					indw = find(lall{w} < point(w), 1, 'last');
					ind{w} = [indw, indw+1];
				elseif indw == 1
					ind{w} = [indw, indw+1];
				else  % indw == end
					ind{w} = [indw-1, indw];
				end
			end
			
			l = cell(1, Axis.count);
			for w = Axis.elems
				l{w} = lall{w}(ind{w});
			end
			[X, Y, Z] = ndgrid(l{:});
			V = this.array(ind{:});
			val = interpn(X, Y, Z, V, point(Axis.x), point(Axis.y), point(Axis.z));
		end
	end		
end
