classdef Grid1d < handle
    % Grid has all the information related to the staggered grid in a single
    % Cartesian axis.  It does not have physical quantities dependent on
    % frequencies, e.g., omega, eps, mu, and PML s-factors.

    properties (SetAccess = immutable)
		axis  % axis of this grid; one of Axis.x, Axis.y, Axis.z
		unit  % instance of PhysUnit
		unitvalue  %  unit value of length
        l  % {locations for primary vertices, locations for dual vertices}
        dl  % {dl at primary vertex, dl at dual vertex} == { diff( l{GK.dual} ), diff( l{GK.prim} ) }
        bc  % [BC at (-) end, BC at (+) end]
        N  % number of grid cells
		L  % length of the grid (domain)
        Npml  % [Npml at (-) end, Npml at (+) end]
		lpml  % [location of PML interface at (-) end, location of PML interface at (+) end]
		Lpml  % [thickness of PML interface at (-) end, thickness of PML interface at (+) end]
		center  % center of grid without PML
	end

	properties (Access = private)
		lghost  % [location of ghost primary vertex (after the last one), location of ghost dual vertex (before the first one)]
		ldual_ext  % location of dual vertex after the last one
	end

	properties (Dependent, SetAccess = immutable)
		lg  % {l_prim with ghost, l_dual with ghost}: l with ghost vertex (beyond boundary)
		lall  % {l_prim with ghost, l_dual with extra vertices}: l with vertices to interpolate fields at corners of simulation domain
		bound  % [lall_prim(1), lall_prim(end)]
	end
	
	properties (SetAccess = private)
		kBloch
	end
        
    methods
        function this = Grid1d(axis, unit, lprim_array, Npml_array, bc)
			chkarg(istypesizeof(axis, 'Axis'), '"axis" should be instance of Axis.');
			this.axis = axis;
			
			chkarg(istypesizeof(unit, 'PhysUnit'), '"unit" should be instance of PhysUnit.');
			this.unit = unit;
			this.unitvalue = unit.value(PhysQ.L);
			
			chkarg(istypesizeof(Npml_array, 'int', [1, Sign.count]), ...
				'"Npml" should be length-%d row vector with integer elements.', Sign.count);
			this.Npml = Npml_array;
			            
			chkarg(istypesizeof(bc, 'BC'), '"bc" should be instance of BC.');
			if bc == BC.p
				this.bc = [bc bc];
			else
				this.bc = [bc BC.m];
			end
			
			chkarg(istypesizeof(lprim_array, 'real', [1 0]), ...
				'"lprim_array" should be row vector with real elements.');
			
			% Set N and L.
			lprim = lprim_array;
			this.N = length(lprim)-1;  % # of grid cells in the axis
			this.L = lprim(end) - lprim(1);
			
			% Set loc and dl.
			ldual = NaN(1, this.N+1);
			ldual(2:(this.N+1)) = (lprim(1:end-1) + lprim(2:end)) / 2;

			if this.bc(Sign.n) == BC.p  % this.bc(Sign.p) == BC.p
				ldual(1) = ldual(end) - (lprim(end)-lprim(1));  % lprim(end) - lprim(1) == ldual(end) - ldual(1)
				this.ldual_ext = ldual(2) + (lprim(end)-lprim(1));  % lprim(end) - lprim(1) = ldual_ext - ldual(2)
			else  % this.bc(Sign.p) == BC.m
				ldual(1) = lprim(1) - (ldual(2)-lprim(1));  % lprim(1) - ldual(1) == ldual(2) - lprim(1)
				this.ldual_ext = lprim(end) + (lprim(end) - ldual(end));  % ldual_ext - lprim(end) = lprim(end) - ldual(end)
			end
			
			this.l = {lprim(1:end-1), ldual(2:end)};
			this.lghost = [lprim(end), ldual(1)];
			this.dl = {diff(ldual), diff(lprim)};  % not {diff(lprim), diff(ldual)}
			
			% Set lpml, Lpml, and center.
			this.lpml = [lprim(1 + this.Npml(Sign.n)), lprim(end - this.Npml(Sign.p))];
			this.Lpml = [this.lpml(Sign.n) - lprim(1), lprim(end) - this.lpml(Sign.p)];
			this.center = mean(this.lpml);
			
			% Initialize kBloch
			this.kBloch = 0;
		end
		
		function lg = get.lg(this)
			lg = {[this.l{GK.prim}, this.lghost(GK.prim)], [this.lghost(GK.dual), this.l{GK.dual}]};
		end
		
		function lall = get.lall(this)
			lall = {this.lg{GK.prim}, [this.lg{GK.dual}, this.ldual_ext]};
		end
		
		function bound = get.bound(this)
			bound = this.lall{GK.prim}([1 end]);
		end
		
		function set_kBloch(this, plane_src)
			chkarg(istypesizeof(plane_src, 'PlaneSrc'), '"plane_src" should be instance of PlaneSrc.');
			this.kBloch = plane_src.kBloch(this.axis);
		end
		
		function truth = contains(this, l)
			% This function can handle "l" as an array.
			chkarg(istypesizeof(l, 'real', [0 1]), '"l" should be column vector with real elements.');
			truth = (l >= this.lall{GK.prim}(1)) & (l <= this.lall{GK.prim}(end));  % &: elementwise AND operator
		end
		
		function bound_plot = bound_plot(this, withpml)
			lplot = this.lplot(GK.prim, withpml);
			bound_plot = lplot([1 end]);
		end
		
		function lplot = lplot(this, gk, withpml)
			chkarg(istypesizeof(gk, 'GK') , '"gk" should be instance of GK');
			chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');
			
			lplot = this.lall{gk};
			if ~withpml
				lplot = lplot(1+this.Npml(Sign.n):end-this.Npml(Sign.p));
			end
		end
	end
end
