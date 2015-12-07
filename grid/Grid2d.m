classdef Grid2d < handle
    % Grid2d has all the information related to the 2D Yee grid.
	% It does not have physical quantities dependent on frequencies, e.g.,
	% omega, eps, mu, and PML s-factors.

    properties (SetAccess = immutable)
		normal_axis
		comp  % [Grid1d_h, Grid1d_v]
		unitvalue  %  unit value of length
    end
        
	properties (Dependent, SetAccess = immutable)
		axis  % [h_axis, v_axis]: for normal_axis == Axis.y, axis = [Axis.z, Axis.x]
		unit  % instance of PhysUnit
        l  % {h_prim, h_dual; v_prim, v_dual}
		lg  %  {h_prim with ghost, h_dual with ghost; v_prim with ghost, v_dual with ghost}
		lall  % {h_prim with ghost, h_dual with extra vertices; v_prim with ghos, v_dual with extra vertices}
		bound  % [hall_prim(1), hall_prim(end); vall_prim(1), vall_prim(end)]
        dl  % {diff(h_dual), diff(h_prim); diff(v_dual), diff(v_prim)}
        bc  % [bc_h, bc_v]
        N  % [Nh, Nv]: # of grid cells in the horizontal and vertical directions
		Ncell  % {[Nh], [Nv]}: to pass as function arguments
		Ntot  % Nh * Nv
		L  % [Lh, Lv]; size of the domain
        Npml  % [Npml_hn, Npml_hp; Npml_vn, Npml_vp]
        lpml  % [lpml_hn, lpml_hp; lpml_vn, lpml_vp]
        Lpml  % [Lpml_hn, Lpml_hp; Lpml_vn, Lpml_vp]
		center  % [center_h, center_v]: center of grid without PML
	end
	
	properties (Dependent, SetAccess = immutable)
		kBloch
	end
        
	methods
        function this = Grid2d(grid3d, normal_axis)
			% Check and store arguments.
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.'); 
			chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
			[h, v] = cycle(normal_axis);
			
			this.normal_axis = normal_axis;
			this.comp = [grid3d.comp(h), grid3d.comp(v)];
		end
		
		function axis = get.axis(this)
			axis = Axis.empty();
			for d = Dir.elems
				axis(d) = this.comp(d).axis;
			end
		end
		
		function unit = get.unit(this)
			unit = this.comp(Dir.h).unit;
		end
		
		function unit = get.unitvalue(this)
			unit = this.comp(Dir.h).unitvalue;
		end
		
		function l = get.l(this)
			l = cell(Dir.count, GT.count);
			for d = Dir.elems
				for g = GT.elems
					l{d, g} = this.comp(d).l{g};
				end
			end
		end
		
		function lg = get.lg(this)
			lg = cell(Dir.count, GT.count);
			for d = Dir.elems
				for g = GT.elems
					lg{d, g} = this.comp(d).lg{g};
				end
			end
		end
			
		function lall = get.lall(this)
			lall = cell(Dir.count, GT.count);
			for d = Dir.elems
				for g = GT.elems
					lall{d, g} = this.comp(d).lall{g};
				end
			end
		end
		
		function bound = get.bound(this)
			bound = NaN(Dir.count, Sign.count);
			for d = Dir.elems
				bound(d,:) = this.comp(d).bound;
			end
		end
		
		function dl = get.dl(this)
			dl = cell(Dir.count, GT.count);
			for d = Dir.elems
				for g = GT.elems
					dl{d, g} = this.comp(d).dl{g};
				end
			end
		end
		
		function bc = get.bc(this)
			bc = BC.empty(0, Dir.count);
			for d = Dir.elems
				bc(d) = this.comp(d).bc;
			end
		end
		
		function N = get.N(this)
			N = NaN(1, Dir.count);
			for d = Dir.elems
				N(d) = this.comp(d).N;
			end
		end
		
		function Ncell = get.Ncell(this)
			Ncell = num2cell(this.N);
		end
		
		function Ntot = get.Ntot(this)
			Ntot = prod(this.N);
		end
		
		function L = get.L(this)
			L = NaN(1, Dir.count);
			for d = Dir.elems
				L(d) = this.comp(d).L;
			end
		end
		
		function Npml = get.Npml(this)
			Npml = NaN(Dir.count, Sign.count);
			for d = Dir.elems
				Npml(d,:) = this.comp(d).Npml;
			end
		end
		
		function lpml = get.lpml(this)
			lpml = NaN(Dir.count, Sign.count);
			for d = Dir.elems
				lpml(d,:) = this.comp(d).lpml;
			end
		end
		
		function Lpml = get.Lpml(this)
			Lpml = NaN(Dir.count, Sign.count);
			for d = Dir.elems
				Lpml(d,:) = this.comp(d).Lpml;
			end
		end
		
		function center = get.center(this)
			center = NaN(1, Dir.count);
			for d = Dir.elems
				center(d) = this.comp(d).center;
			end
		end
		
		function kBloch = get.kBloch(this)
			kBloch = NaN(1, Dir.count);
			for d = Dir.elems
				kBloch(d) = this.comp(d).kBloch;
			end
		end
				
		function truth = contains(this, h, v)
			chkarg(istypeof(h, 'real'), '"h" should be array with real elements.');
			chkarg(istypeof(v, 'real'), '"v" should be array with real elements.');
			chkarg(isequal(size(h), size(v)), '"h" and "v" should have same size.');
			
			truth = true(size(h));
			loc = {h, v};
			for d = Dir.elems
				truth = truth & this.comp(d).contains(loc{d});  % &: elementwise AND operator
			end
		end
		
		function bound_plot = bound_plot(this, withpml)
			bound_plot = NaN(Dir.count, Sign.count);
			for d = Dir.elems
				bound_plot(d,:) = this.comp(d).bound_plot(withpml);
			end
		end
		
		function lplot_cell = lplot(this, g, withinterp, withpml)
			chkarg(istypesizeof(g, 'GT') || istypesizeof(g, 'GT', [1 Dir.count]), '"g" should be instance of GT');
			if length(g) == 1
				g = g(ones(1, Dir.count));
			end
			chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
			chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');
			
			lplot_cell = cell(1, Dir.count);
			for d = Dir.elems
				lplot_cell{d} = this.comp(d).lplot(g(d), withinterp, withpml);
			end
		end
		
		function lvoxelbound_cell = lvoxelbound(this, g, withpml)
			chkarg(istypesizeof(g, 'GT') || istypesizeof(g, 'GT', [1 Dir.count]), '"g" should be instance of GT');
			if length(g) == 1
				g = g(ones(1, Dir.count));
			end
			chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');

			lvoxelbound_cell = cell(1, Dir.count);
			for d = Dir.elems
				lvoxelbound_cell{d} = this.comp(d).lvoxelbound(g(d), withpml);
			end
		end
	end
end
