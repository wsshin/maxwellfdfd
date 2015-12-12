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
		
		function l_cell = lall(this)
			l_cell = this.grid2d.lall(Dir.elems + Dir.count*subsindex(this.gt_array));
		end
		
		function l_cell = lg(this)
			l_cell = this.grid2d.lg(Dir.elems + Dir.count*subsindex(this.gt_array));
		end
		
		function l_cell = l(this)
			l_cell = this.grid2d.l(Dir.elems + Dir.count*subsindex(this.gt_array));
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
		
		function val = value(this, h, v)
			loc = {h, v};
			lavail = this.l;  % available locations
			num_len1 = 0;
			for d = Dir.elems
				coord_d = loc{d};
				chkarg(isempty(coord_d) || isvector(coord_d) && istypeof(coord_d, 'real'), ...
				'"%s" should be empty or real vector.', char(d));
				if isempty(coord_d)
					coord_d = lavail{d};
				else
					chkarg(this.grid2d.comp(d).contains(coord_d), '"%s" should be inside grid.', char(d));
					coord_d = sort(coord_d);
				end
				loc{d} = coord_d;
				if length(coord_d) == 1
					num_len1 = num_len1 + 1;
				end
			end
			chkarg(num_len1 >= 1, 'Points should be along Cartesian direction.');
			
			lall = this.lall;
			ind = cell(1, Dir.count);
			
			for d = Dir.elems
				loc_min = loc{d}(1);
				indw_min = ismembc2(loc_min, lall{d});
				if indw_min == 0
					indw_min = find(lall{d} < loc_min, 1, 'last');
				end
				
				loc_max = loc{d}(end);
				indw_max = ismembc2(loc_max, lall{d});
				if indw_max == 0
					indw_max = find(lall{d} > loc_max, 1, 'first');
				end
				
				ind{d} = indw_min:indw_max;
			end
			
			isindlen1 = [length(ind{Dir.h}), length(ind{Dir.v})] == 1;
			if all(isindlen1)
				val = this.array(ind{:});
			else
				for d = Dir.elems
					if isindlen1(d) % ndgrid() needs two data points to interpolate with
						indw = ind{d};
						if indw >= 2
							ind{d} = [indw - 1, indw];
						else
							ind{d} = [indw, indw + 1];
						end
					end
				end
				
				l = cell(1, Dir.count);
				for d = Dir.elems
					l{d} = lall{d}(ind{d}(1):ind{d}(end));
				end
				[H, V] = ndgrid(l{:});
				C = this.array(ind{:});
				[Hi, Vi] = ndgrid(loc{:});
				val = interpn(H, V, C, Hi, Vi);
			end			
		end
	end		
end
