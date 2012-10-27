% Note that wvlen can be different from Oscillation.in_L0() because wvlen here
% should be the wavelength inside the medium, not inside a vacuum.

classdef PlaneSrc < Source
	% PlaneSrc is a class representing a constant electric dipole distribution
	% over the entire plane.  PlaneSrc is to generate a plane wave in a
	% homogeneous medium, and it supports the Bloch boundary condition.
	
	properties (SetAccess = immutable)
		normal_axis  % plane normal axis: one of Axis.x, Axis.y, Axis.z
		intercept  % intercept between plane and normal axis
		phi  % angle of polarization with respect to first Cartesian direction in this plane
		theta  % emission angle in the plane normal to polarization: (+) toward normal x polarization, (-) toward -normal x polarization
		kBloch  % Bloch k-vector [kx, ky, kz]
		K  % surface current density
	end
	
	methods
		function this = PlaneSrc(normal_axis, intercept, polarization, K, theta, wvlen)
			chkarg(istypesizeof(normal_axis, 'Axis'), ...
				'"normal_axis" should be instance of Axis.');
			chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real.');
			
			chkarg(istypesizeof(polarization, 'Axis') || istypesizeof(polarization, 'real'), ...
				'"polarization" should be instance of Axis or angle in radian.');
			if istypesizeof(polarization, 'Axis')
				[p, q] = cycle(normal_axis);
				if polarization == p
					polarization = 0;
				else
					assert(polarization == q);
					polarization = pi/2;
				end
			end
			
			if nargin < 4  % no K
				K = 1.0;
			end
			chkarg(istypesizeof(K, 'complex'), '"K" should be complex.');

			if nargin < 6  % no wvlen
				assert(nargin < 5);  % no theta
				theta = 0;
				wvlen = inf;
			end
			chkarg(istypesizeof(theta, 'real'), '"theta" should be polar angle in radian.');
			chkarg(istypesizeof(wvlen, 'real') && wvlen > 0, '"wvlen" should be positive.');
						
			l = cell(Axis.count, GK.count);
			l{normal_axis, GK.dual} = intercept;
			this = this@Source(l);
			
			this.normal_axis = normal_axis;
			this.intercept = intercept;
			this.phi = polarization;
			this.theta = theta;
			this.K = K;

			this.kBloch = [0 0 0];
			if this.theta ~= 0
				[p, q] = cycle(this.normal_axis);
				kt = (2*pi/wvlen) * sin(this.theta);  % k tangential to plane
				kp = kt * cos(this.phi + pi/2);
				kq = kt * sin(this.phi + pi/2);
				this.kBloch(p) = kp;
				this.kBloch(q) = kq;
			end
		end
		
		function [index_cell, Jw_patch] = generate_kernel(this, w_axis, grid3d)
			grid3d.set_kBloch(this);
			index_cell = cell(1, Axis.count);
			if w_axis == this.normal_axis
				Jw_patch = [];
			else
				n = this.normal_axis;
				w = w_axis;
				v = setdiff(Axis.elems, [n, w]);
				
				g = GK.dual;
				ind_n = ismembc2(this.intercept, grid3d.l{n,g});
				if ind_n == 0
					[~, ind_n] = min(abs(grid3d.l{n,g} - this.intercept));
					warning('FDS:srcAssign', ...
						['%s grid in %s-axis of "grid3d" does not have location %e of this %s; ', ...
						'closest grid vertex at %e will be taken instead.'], ...
						char(g), char(n), this.intercept, class(this), grid3d.l{n,g}(ind_n));
				end

				dn = grid3d.dl{n,g}(ind_n);
				J = this.K / dn;  % for t normal to n and K, K*dt = (current through dn*dt) = J * (dn*dt)
				
				if w == cycle(this.normal_axis)
					Jw = J * cos(this.phi);
				else
					Jw = J * sin(this.phi);
				end
				
				if Jw == 0
					Jw_patch = [];
				else
					%  Set Jw_patch.
					lw = grid3d.l{w, GK.prim};
					lv = grid3d.l{v, GK.dual};
					kw = this.kBloch(w);
					kv = this.kBloch(v);

% 					if w < v
% 						Jw_patch = Jw .* (exp(-1i * (kw .* lw)).' * exp(-1i * (kv .* lv)));
% 					else
% 						assert(w > v);
% 						Jw_patch = Jw .* (exp(-1i * (kv .* lv)).' * exp(-1i * (kw .* lw)));
% 					end

					Jw_patch = Jw .* (exp(-1i * (kw .* lw)).' * exp(-1i * (kv .* lv)));
					Jw_patch = ipermute(Jw_patch, int([w v n]));

					% Set index_cell.
					index_cell{n} = ind_n;
					Nw = grid3d.N(w);
					Nv = grid3d.N(v);
					index_cell{w} = 1:Nw;
					index_cell{v} = 1:Nv;
				end
			end
		end		
	end
end

