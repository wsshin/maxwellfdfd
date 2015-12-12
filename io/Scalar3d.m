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
		
		function l_cell = lall(this)
			l_cell = this.grid3d.lall(Axis.elems + Axis.count*subsindex(this.gt_array));
		end
		
		function l_cell = lg(this)
			l_cell = this.grid3d.lg(Axis.elems + Axis.count*subsindex(this.gt_array));
		end
		
		function l_cell = l(this)
			l_cell = this.grid3d.l(Axis.elems + Axis.count*subsindex(this.gt_array));
		end

		function [array, l_cell] = data_expanded(this)
			l_cell = this.lall;
			array = this.array;
		end
		
		function [array, l_cell] = data_ghost_expanded(this)
			l_cell = this.lg;
			ind = cell(1, Axis.count);
			
			for w = Axis.elems
				ind{w} = 1:(this.grid3d.N(w)+1);
			end
			array = this.array(ind{:});
		end
		
		function [array, l_cell] = data_original(this)
			l_cell = this.l;
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
		
		function [val, loc] = value(this, x, y, z)
			% Return the values of this Scalar3d along a line.  The line should
			% be along one of the x, y, z directions, so at least two of the
			% given x, y, z should be scalar.
			loc = {x, y, z};
			lavail = this.l;  % available locations
			num_len1 = 0;  % at least two of x, y, z should be scalar
			for w = Axis.elems
				coord_w = loc{w};
				chkarg(isempty(coord_w) || isvector(coord_w) && istypeof(coord_w, 'real'), ...
				'"%s" should be empty or real vector.', char(w));
				if isempty(coord_w)
					coord_w = lavail{w};  % if empty, use all available locations
				else
					chkarg(this.grid3d.comp(w).contains(coord_w), '"%s" should be inside grid.', char(w));
					coord_w = sort(coord_w);
				end
				loc{w} = coord_w;
				if length(coord_w) == 1
					num_len1 = num_len1 + 1;
				end
			end
			chkarg(num_len1 >= 2, 'Points should be along Cartesian direction.');
			
			lall = this.lall;
			ind = cell(1, Axis.count);
			
			for w = Axis.elems
				loc_min = loc{w}(1);
				indw_min = ismembc2(loc_min, lall{w});
				if indw_min == 0
					indw_min = find(lall{w} < loc_min, 1, 'last');
				end
				
				loc_max = loc{w}(end);
				indw_max = ismembc2(loc_max, lall{w});
				if indw_max == 0
					indw_max = find(lall{w} > loc_max, 1, 'first');
				end
				
				ind{w} = indw_min:indw_max;
			end
			
			isindlen1 = [length(ind{Axis.x}), length(ind{Axis.y}), length(ind{Axis.z})] == 1;
			if all(isindlen1)
				val = this.array(ind{:});
			else
				for w = Axis.elems
					if isindlen1(w) % ndgrid() needs two data points to interpolate with
						indw = ind{w};
						if indw >= 2
							ind{w} = [indw - 1, indw];
						else
							ind{w} = [indw, indw + 1];
						end
					end
				end
				
				l = cell(1, Axis.count);
				for w = Axis.elems
					l{w} = lall{w}(ind{w}(1):ind{w}(end));
				end
				[X, Y, Z] = ndgrid(l{:});
				V = this.array(ind{:});
				[Xi, Yi, Zi] = ndgrid(loc{:});
				val = interpn(X, Y, Z, V, Xi, Yi, Zi);
			end
		end
	end		
end
