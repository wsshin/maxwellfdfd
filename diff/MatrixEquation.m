% Oscillation is a class containing the information on the operation frequency.
classdef MatrixEquation
    properties (SetAccess = immutable, GetAccess = private)
	end

    properties (SetAccess = immutable)
		Ce  % matrix of curl for E-field
		Cm  % matrix of curl for H-field
		EPS  % matrix of subpixel-smoothed permittivity tensor
		MU  % matrix of subpixel-smoothed permeability tensor
		j  % column vector of electric current source
		m  % column vector of magnetic current source
		omega  % angular frequency
		pecm  % column vector of PEC mask
		r  % reordering array
		s  % inverse reordering array
		ft  % field type of equation
    end
	
	methods
        function this = MatrixEquation(eqtype, pml, omega, eps_array, mu_array, s_factor_cell, J_cell, M_cell, grid3d)
			chkarg(istypesizeof(eqtype, 'EquationType'), '"eqtype" should be instance of EquationType.');
			chkarg(istypesizeof(pml, 'PML'), '"pml" should be instance of PML.');
			chkarg(istypesizeof(omega, 'complex'), '"omega" should be complex.');
			chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

			chkarg(istypesizeof(eps_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
				'"eps_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});
			chkarg(istypesizeof(mu_array, 'complex', [Axis.count Axis.count grid3d.N]), ...
				'"mu_array" should be %d-by-%d-by-%d-by-%d-by-%d array with complex elements.', Axis.count, Axis.count, grid3d.Ncell{:});
			chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
				'"mu_cell" should be %d-by-%d cell array whose each element is row vector with complex elements', ...
				Axis.count, GT.count);
			chkarg(istypesizeof(J_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
				'"J_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
				Axis.count, grid3d.Ncell{:});
			chkarg(istypesizeof(M_cell, 'complexcell', [1 Axis.count], grid3d.N), ...
				'"M_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
				Axis.count, grid3d.Ncell{:});

			% Reorder the indices of the elements of matrices and vectors to reduce the bandwidth of A.
			this.r = reordering_indices(Axis.count, grid3d.N);
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
				gm = alter(ge);
				se = s_factor_cell(:, ge).';
				sm = s_factor_cell(:, gm).';
				
				for w = Axis.elems
					for v = Axis.elems
						sdim = {1, 1, 1};
						sdim{v} = grid3d.N(v);
					
						if v ~= w  % left-multiply area elements to eps and mu tensors, so take rows of eps and mu
							eps_array(w,:,:,:,:) = bsxfun(@times, eps_array(w,:,:,:,:), reshape(se(v), 1, 1, sdim{:}));
							mu_array(w,:,:,:,:) = bsxfun(@times, mu_array(w,:,:,:,:), reshape(sm(v), 1, 1, sdim{:}));
						else  % v == w; right-multiply inverse of length elements to eps and mu tensor, so take columns of eps and mu
							eps_array(:,w,:,:,:) = bsxfun(@rdivide, eps_array(:,w,:,:,:), reshape(sm(v), 1, 1, sdim{:}));
							mu_array(:,w,:,:,:) = bsxfun(@rdivide, mu_array(:,w,:,:,:), reshape(se(v), 1, 1, sdim{:}));
						end
					end
				end
			end

			[this.EPS, this.MU] = create_EPS_MU(eps_array, mu_array, ge, dl_factor_cell, grid3d);

			this.j = [J_cell{Axis.x}(:) ; J_cell{Axis.y}(:) ; J_cell{Axis.z}(:)];
			this.m = [M_cell{Axis.x}(:) ; M_cell{Axis.y}(:) ; M_cell{Axis.z}(:)];
			
			this.omega = omega;

			% Mask elements corresponding to PEC.
			ind_pec = isinf(max(abs(this.EPS)));
			
			% We know the E-fields at ind_pec locations are all zero, so we
			% could a new set of linear equations witout the corresponding
			% columns in the matrix, which will give an overdetermined set of
			% linear equations.  We can recover a square matrix by eliminating
			% the corresponding rows as well, because we should be able to solve
			% the equations without the rows anyway.  Now, instead of reducing
			% the size of the matrix, we just leave only the diagonal entries in
			% the corresponding columns and rows.  This ensures the
			% corresponding E-field to be zero when the corresponding J is zero.
			this.EPS(:, ind_pec) = 0;
			this.EPS(ind_pec, :) = 0;
			this.EPS = this.EPS + create_spdiag(double(ind_pec));
			
			this.pecm = double(~ind_pec);  % PEC mask

			test_vector_identity = false;
			if test_vector_identity
				[Dive, Divm] = create_divs(ge, dl_factor_cell, grid3d);
				fprintf('norm(Dive * Cm, 1) = %e\n', norm(Dive * this.Cm, 1));
				fprintf('norm(Divm * Ce, 1) = %e\n', norm(Divm * this.Ce, 1));
			end
			
			this.ft = eqtype.f;
		end
		
        function [A, b] = matrix_op(this)
			if this.ft == FT.e
% 				INV_MU = create_spdiag(1./this.mu);  % create INV_MU instead of inverting MU; "MU \ Mat" complains about singularity when mu has Inf
% 				EPS = create_spdiag(this.eps);
				PECM = create_spdiag(this.pecm);	

				A = PECM * (this.Cm * (this.MU \ this.Ce)) * PECM - this.omega^2 * this.EPS;
				b = -1i * this.omega * this.j - this.Cm * (this.MU \ this.m);
            
				A = A(this.r, this.r);
				b = b(this.r);
			elseif this.ft == FT.h
% 				INV_EPS = create_spdiag(1./this.eps);  % create INV_EPS instead of inverting EPS; "EPS \ Mat" complains about singularity when eps has Inf
% 				MU = create_spdiag(this.mu);

				A = (this.Ce * (this.EPS \ this.Cm)) - this.omega^2 * this.MU;
				b = -1i * this.omega * this.m + this.Ce * (this.EPS \ this.j);
            
				A = A(this.r, this.r);
				b = b(this.r);
			else  % use both E and H
				assert(isequal(this.ft, FT.elems));
				PECM = create_spdiag(this.pecm);	
				
				A = [-1i * this.omega * this.EPS, PECM * this.Cm; this.Ce * PECM, 1i * this.omega * this.MU];
				b = [this.j; -this.m];
            
				N = length(this.r);
				rr = [this.r; (this.r + N)];
				A = A(rr, rr);
				b = b(rr);
			end
		end
		
		function [h_Op, b, h_GfromF] = matrixfree_op(this)
			if this.ft == FT.e
				% A = PECM * (this.Cm * (this.MU \ this.Ce)) * PECM - this.omega^2 * this.EPS;
				h_A = @(e) this.pecm .* (this.Cm * (this.MU \ (this.Ce * (this.pecm .* e)))) - this.omega^2 * (this.EPS * e);
				h_Atr = @(e) this.pecm .* (this.Ce.' * (this.MU.' \ (this.Cm.' * (this.pecm .* e)))) - this.omega^2 * (this.EPS.' * e);

				b = -1i * this.omega * this.j - this.Cm * (this.MU \ this.m);

				h_GfromF_temp = @(e) (this.MU \ (this.Ce * e + this.m)) ./ (-1i*this.omega);
			else  % this.ft == FT.h
				% A = (this.Ce * (this.EPS \ this.Cm)) - this.omega^2 * this.MU;
				h_A = @(h) this.Ce * (this.EPS \ (this.Cm * h)) - this.omega^2 * (this.MU * h);
				h_Atr = @(h) this.Cm.' * (this.EPS.' \ (this.Ce.' * h)) - this.omega^2 * (this.MU.' * h);

				b = -1i * this.omega * this.m + this.Ce * (this.EPS \ this.j);

				h_GfromF_temp = @(h) (this.EPS \ (this.Cm * h - this.j)) ./ (1i*this.omega);
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

