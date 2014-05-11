classdef Scalar2d
    % SCALAR2D represents a 2D arary of scalar values.  It can represents a single
    % directional component of a 3D vector field (e.g. Ex of E).
    
    properties (SetAccess = immutable)
		% numerical properties
        array  % 2D array
        grid2d  % instance fo Grid2d
		gt_array  % [h_gridtype v_gridtype]: array of instances of GT
        osc  % instance of Oscillation
		intercept  % intercept of the plane of this 2D data in the normal axis

		% extra properties
		physQcell  % {PhysQ, int; PhysQ, int; ...}: representation of physical quanity
		unitvalue  % unit value of physical quantity
		name  % name of this quantity: "E_x", "H_y", etc. in LaTeX equation format
	end
	            
    methods
        function this = Scalar2d(array, grid2d, gt_array, osc, physQcell, name, intercept)
			chkarg(istypesizeof(grid2d, 'Grid2d'), '"grid2d" should be instance of Grid2d.');
			this.grid2d = grid2d;
			
			N = grid2d.N + int(gt_array);
			chkarg(istypesizeof(array, 'complex', N), ...
				'"array" should be %d-by-%d array with complex elements.', N(Dir.h), N(Dir.v));
			this.array = array;
			
			chkarg(istypesizeof(gt_array, 'GT', [1 Dir.count]), ...
				'"gt_array" should be length-%d row vector with GT as elements.', Dir.count);
			this.gt_array = gt_array;
			
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
			
			if nargin < 7  % no intercept
				intercept = NaN;
			end
			this.intercept = intercept;
		end
		
		function l_cell = lplot(this, withinterp, withpml)
			% Return the locations where data are evaluated for plotting.  If
			% the data do not include the boundaries of the simulation domain (or
			% the PML interfaces for "withpml == false"), the boundary points
			% are added.
			l_cell = this.grid2d.lplot(this.gt_array, withinterp, withpml);
		end
		
		function l_cell = lvoxelbound(this, withpml)
			% Return the locations of boundaries of voxels drawn.  For data at
			% primary grid points, the voxel centers are in the simulation
			% domain including the boundary.  For data at dual grid points, the
			% voxel centers are in the simulation domain excluding the boundary.
			l_cell = this.grid2d.lvoxelbound(this.gt_array, withpml);
		end
		
		function [array, l_cell] = data_expanded(this)
			l_cell = this.grid2d.lall(Dir.elems + Dir.count*subsindex(this.gt_array));
			array = this.array;
		end
		
		function [array, l_cell] = data_original(this)
			l_cell = this.grid2d.l(Dir.elems + Dir.count*subsindex(this.gt_array));
			ind = cell(1, Dir.count);
			
			for d = Dir.elems
				if this.gt_array(d) == GT.prim
					ind{d} = 1:this.grid2d.N(d);
				else  % this.gt_array(d) == GT.dual
					ind{d} = 2:(this.grid2d.N(d)+1);
				end
			end
			array = this.array(ind{:});
		end
		
		function [C, X, Y] = data_for_pcolor(this, withinterp, withpml)
			chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
			chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');
			
			lall = this.grid2d.lall(Dir.elems + Dir.count*subsindex(this.gt_array));
			[Xh, Yv] = ndgrid(lall{:});

			lplot = this.lplot(withinterp, withpml);
			[XIh, YIv] = ndgrid(lplot{:});
			C = interpn(Xh, Yv, this.array, XIh, YIv);
			C = permute(C, int([Dir.v, Dir.h]));  % to be compatible with pcolor()

			if ~withinterp
				% Pad the end boundaries.  The padded values are not drawn, so
				% they don't have to be accurate.
				C(end+1, :) = C(end, :);
				C(:, end+1) = C(:, end);
			end
			
			if nargout >= 2  % X, Y
				if withinterp
					l = this.lplot(withinterp, withpml);
				else
					l = this.lvoxelbound(withpml);
				end
				[X, Y] = meshgrid(l{:});
			end
		end
			
		function val = value(this, point)
			chkarg(istypesizeof(point, 'real', [0, Dir.count]), ...
				'"point" should be 2D array with %d columns.', Dir.count);
			chkarg(this.grid2d.contains(point), '"point" should be inside grid.');
			
			lall = this.grid2d.lall(Dir.elems + Dir.count*subsindex(this.gt_array));
			ind = cell(1, Dir.count);
			for d = Dir.elems
				indd = ismembc2(point(d), lall{d});
				if indd == 0
					indd = find(lall{d} < point(d), 1, 'last');
					ind{d} = [indd, indd+1];
				elseif indd == 1
					ind{d} = [indd, indd+1];
				else  % indd == end
					ind{d} = [indd-1, indd];
				end
			end
			
			l = cell(1, Dir.count);
			for d = Dir.elems
				l{d} = lall{d}(ind{d});
			end
			[X, Y] = ndgrid(l{:});
			C = this.array(ind{:});
			val = interpn(X, Y, C, point(Dir.h), point(Dir.v));
		end		
	end		
end
