% Oscillation is a class containing the information on the operation frequency.
classdef MatrixEquation
    properties (SetAccess = immutable, GetAccess = private)
	end

    properties (SetAccess = immutable)
		Ce  % matrix of curl for E-field
		Cm  % matrix of curl for H-field
		eps  % column vector of permittivity
		mu  % column vector of permeability
		j  % column vector of electric current source
		m  % column vector of magnetic current source
		omega  % angular frequency
		pm  % column vector of PEC mask
		r  % reordering array
		s  % inverse reordering array
		ft  % field type of equation
    end
	
	methods
        function this = MatrixEquation(eqtype, pml, omega, eps_cell, mu_cell, s_factor_cell, J_cell, M_cell, grid3d)
			chkarg(istypesizeof(eqtype, 'EquationType'), '"eqtype" should be instance of EquationType.');
			chkarg(istypesizeof(pml, 'PML'), '"pml" should be instance of PML.');
			chkarg(istypesizeof(omega, 'complex'), '"omega" should be complex.');
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

			N = grid3d.N;
			chkarg(istypesizeof(eps_cell, 'complexcell', [1 Axis.count], N), ...
				'"eps_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
				Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));
			chkarg(istypesizeof(mu_cell, 'complexcell', [1 Axis.count], N), ...
				'"mu_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
				Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));
			chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
				'"mu_cell" should be %d-by-%d cell array whose each element is row vector with complex elements', ...
				Axis.count, GT.count);
			chkarg(istypesizeof(J_cell, 'complexcell', [1 Axis.count], N), ...
				'"J_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
				Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));
			chkarg(istypesizeof(M_cell, 'complexcell', [1 Axis.count], N), ...
				'"M_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
				Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));

			% Reorder the indices of the elements of matrices and vectors to reduce the bandwidth of A.
			this.r = reordering_indices(Axis.count, N);
			[~, this.s] = sort(this.r);  % s is inverse reordering

			% Construct curls
			dl_factor_cell = [];
			if pml == PML.sc
				dl_factor_cell = s_factor_cell;
			end

			ge = eqtype.ge;
			[this.Ce, this.Cm] = create_curls(ge, dl_factor_cell, grid3d);

			% Construct material parameters
			if pml == PML.u
				[smx, smy, smz] = ndgrid(s_factor_cell{Axis.x, alter(ge)}, s_factor_cell{Axis.y, alter(ge)}, s_factor_cell{Axis.z, alter(ge)});
				[sex, sey, sez] = ndgrid(s_factor_cell{Axis.x, ge}, s_factor_cell{Axis.y, ge}, s_factor_cell{Axis.z, ge});
				sm = {smx, smy, smz};
				se = {sex, sey, sez};
				mu_cell = mult_vec(mu_cell, sm([Axis.y Axis.z Axis.x]));
				mu_cell = mult_vec(mu_cell, sm([Axis.z Axis.x Axis.y]));
				mu_cell = div_vec(mu_cell, se([Axis.x Axis.y Axis.z]));

				eps_cell = mult_vec(eps_cell, se([Axis.y Axis.z Axis.x]));
				eps_cell = mult_vec(eps_cell, se([Axis.z Axis.x Axis.y]));
				eps_cell = div_vec(eps_cell, sm([Axis.x Axis.y Axis.z]));
			end

			this.mu = [mu_cell{Axis.x}(:) ; mu_cell{Axis.y}(:) ; mu_cell{Axis.z}(:)];
			this.eps = [eps_cell{Axis.x}(:) ; eps_cell{Axis.y}(:) ; eps_cell{Axis.z}(:)];

			this.j = [J_cell{Axis.x}(:) ; J_cell{Axis.y}(:) ; J_cell{Axis.z}(:)];
			this.m = [M_cell{Axis.x}(:) ; M_cell{Axis.y}(:) ; M_cell{Axis.z}(:)];
			
			this.omega = omega;
			
			this.ft = eqtype.f;

			% Mask elements corresponding to PEC.
			if this.ft == FT.e
				ind_pec = isinf(abs(this.eps));
				this.eps(ind_pec) = 1;
				this.pm = ones(size(ind_pec));  % PEC mask
				this.pm(ind_pec) = 0;
			end

			test_vector_identity = false;
			if test_vector_identity
				[Dive, Divm] = create_divs(ge, dl_factor_cell, grid3d);
				fprintf('norm(Dive * Cm, 1) = %e\n', norm(Dive * this.Cm, 1));
				fprintf('norm(Divm * Ce, 1) = %e\n', norm(Divm * this.Ce, 1));
			end
		end
		
        function [A, b] = matrix_op(this)
			if this.ft == FT.e
				INV_MU = create_spdiag(1./this.mu);  % create INV_MU instead of inverting MU; "MU \ Mat" complains about singularity when mu has Inf
				EPS = create_spdiag(this.eps);
				PM = create_spdiag(this.pm);	

				A = PM * (this.Cm * INV_MU * this.Ce) * PM - this.omega^2 * EPS;
				b = -1i*this.omega*this.j - this.Cm*(this.m./this.mu);
            
				A = A(this.r, this.r);
				b = b(this.r);
			elseif this.ft == FT.h
				INV_EPS = create_spdiag(1./this.eps);  % create INV_EPS instead of inverting EPS; "EPS \ Mat" complains about singularity when eps has Inf
				MU = create_spdiag(this.mu);

				A = (this.Ce * INV_EPS * this.Cm) - this.omega^2 * MU;
				b = -1i*this.omega*this.m + this.Ce*(this.j./this.eps);
            
				A = A(this.r, this.r);
				b = b(this.r);
			else  % use both E and H
				assert(isequal(this.ft, FT.elems));
				EPS = create_spdiag(this.eps);
				MU = create_spdiag(this.mu);
				
				A = [-1i * this.omega * EPS, this.Cm; this.Ce, 1i * this.omega * MU];
				b = [this.j; -this.m];
            
				N = length(this.r);
				rr = [this.r; (this.r + N)];
				A = A(rr, rr);
				b = b(rr);
			end
		end
		
		function [h_Op, b, h_GfromF] = matrixfree_op(this)
			if this.ft == FT.e
				h_A = @(e) this.pm .* (this.Cm * ((this.Ce * (this.pm .* e)) ./ this.mu)) - this.omega^2 * (this.eps .* e);
				h_Atr = @(e) this.pm .* (this.Ce.' * ((this.Cm.' * (this.pm .* e)) ./ this.mu)) - this.omega^2 * (this.eps .* e);

				b = -1i*this.omega*this.j - this.Cm*(this.m./this.mu);

				h_GfromF_temp = @(e) (this.Ce * e + this.m) ./ (-1i*this.omega*this.mu);
			else  % this.ft == FT.h
				h_A = @(h) this.Ce * ((this.Cm * h) ./ this.eps) - this.omega^2 * (this.mu .* h);
				h_Atr = @(h) this.Cm.' * ((this.Ce.' * h) ./ this.eps) - this.omega^2 * (this.mu .* h);

				b = -1i*this.omega*this.m + this.Ce*(this.j./this.eps);

				h_GfromF_temp = @(h) (this.Cm * h - this.j) ./ (1i*this.omega*this.eps);
			end

			function y = Op(x, transp_flag)
				if strcmp(transp_flag,'transp')
				   y = h_Atr(x(this.s));
				elseif strcmp(transp_flag,'notransp')
				   y = h_A(x(this.s));
				end
			   
				y = y(this.r);
			end
			
			function y = GfromF(x)
				y = h_GfromF_temp(x(this.s));
				y = y(this.r);
			end

			h_Op = @Op;
			h_GfromF = @GfromF;
		end
	end
	
end

