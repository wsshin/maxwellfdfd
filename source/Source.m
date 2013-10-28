classdef Source < handle & matlab.mixin.Heterogeneous
	% Source is a superclass for all electric current sources J.
	
	properties (SetAccess = immutable)
		lgrid  % {x_array, y_array, z_array}: locations (intercepts with normal axis) of grid planes of grid type "this.gt"
		laltgrid  % {x_array, y_array, z_array}: locations (intercepts with normal axis) of grid planes of grid type "alter(this.gt)"
		shape  % shape of source; used to draw source
		forceprim  % true or false: flag to control the behavior of get.l()s
	end
	
	properties (SetAccess = private)
		gt  % GT.prim or GT.dual: type of grid line along which dipole sources of this source are defined
	end

	properties (Dependent, SetAccess = immutable)
		l
	end

	methods (Abstract = true)
		[index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
	end
	
	methods
		function this = Source(lgrid_cell, laltgrid_cell, shape, forceprim)
			chkarg(istypesizeof(lgrid_cell, 'realcell', [1 Axis.count], [1 0]), ...
				'"lgrid_cell" should be length-%d row cell array whose each element is row vector with real elements.', Axis.count);
			chkarg(istypesizeof(laltgrid_cell, 'realcell', [1 Axis.count], [1 0]), ...
				'"laltgrid_cell" should be length-%d row cell array whose each element is row vector with real elements.', Axis.count);
			chkarg(istypesizeof(shape, 'Shape'), '"shape" should be instance of Shape.');
			
			if nargin < 4  % no forceprim
				forceprim = false;
			end
			chkarg(istypesizeof(forceprim, 'logical'), '"forceprim" should be logical.');
			
			this.lgrid = lgrid_cell;
			this.laltgrid = laltgrid_cell;
			this.shape = shape;
			this.gt = GT.empty();
			this.forceprim = forceprim;
		end
		
		function set_gridtype(this, gt)
			% If this is SRCJ and E-field grid is primary, then gt = GT.prim.
			% If this is SRCJ and E-field grid is dual, then gt = GT.dual.
			% If this is SRCM and E-field grid is primary, then gt = GT.dual.
			% If this is SRCM and E-field grid is dual, then gt = GT.prim.
			chkarg(istypesizeof(gt, 'GT'), '"gt" should be instance of GT.');
			this.gt = gt;
		end
		
		function l = get.l(this)
			l = cell(Axis.count, GT.count);
			if this.forceprim
				l(:, GT.prim) = this.lgrid.';
				l(:, GT.dual) = this.laltgrid.';
			else
				l(:, this.gt) = this.lgrid.';
				l(:, alter(this.gt)) = this.laltgrid.';
			end
		end
		
		function [index_cell, JMw_patch] = generate(this, w_axis, grid3d)
			chkarg(istypesizeof(w_axis, 'Axis'), '"w_axis" should be instance of Axis.');
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			
			try
				[index_cell, JMw_patch] = this.generate_kernel(w_axis, grid3d);  % Cw_patch: current source
			catch err
				exception = MException('Maxwell:srcAssign', 'Source assignment failed.');
				throw(addCause(exception, err));
			end
		end
	end
end

