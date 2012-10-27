classdef Source < handle & matlab.mixin.Heterogeneous
	% Source is a superclass for all electric current sources J.
	
	properties (SetAccess = immutable)
		l  % {x_prim_array, x_dual_array; y_prim_array, y_dual_array; z_prim_array, z_dual_array}; used to generate grid
	end
	
	methods (Abstract = true)
		[index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
	end
	
	methods
		function this = Source(l_cell)
			chkarg(istypesizeof(l_cell, 'realcell', [Axis.count, GK.count], [1 0]), ...
				'"l_cell" should be %d-by-%d cell array whose each element is row vector with real elements.', Axis.count, GK.count);
			this.l = l_cell;
		end
		
		function [index_cell, Jw_patch] = generate(this, w_axis, grid3d)
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
			chkarg(istypesizeof(w_axis, 'Axis'), '"w_axis" should be instance of Axis.');
			
			try
				[index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d);
			catch err
				exception = MException('FDS:srcAssign', 'Source assignment failed: %s', err.message);
				throw(exception);  % replace with rethrow(err) to reveal err
% 				rethrow(err);
			end
		end
	end
end

