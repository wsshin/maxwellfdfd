function Ds_cell = create_Ds(s, ge, dl_factor_cell, gridnd)
% Creates the forward or backward derivative matrix.

chkarg(istypesizeof(s, 'Sign'), '"s" should be instance of Sign.');  % Ds = Df (or Db) when s == Sign.p (or Sign.n)
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for the E-field
chkarg(istypesizeof(gridnd, 'Grid2d') || istypesizeof(gridnd, 'Grid3d'), '"gridnd" should be instance of Grid2d or Grid3d.');

v = Axis.x;
if istypesizeof(gridnd, 'Grid2d')
	v = Dir.h;
end
chkarg(isempty(dl_factor_cell) || istypesizeof(dl_factor_cell, 'complexcell', [v.count GT.count], [1 0]), ...
	'"dl_factor_cell" should be empty, or %d-by-%d cell array whose each element is row vector with real elements.', v.count, GT.count);

%% Get the shape.
N = gridnd.N;

%% Get the relevant derivative matrices.
g = GT.elems(s);  % curl, divergence, gradient all uses dl_dual for forward difference and dl_prim for backward difference
bc = gridnd.bc;

% Basic setup of Df and Db.  Here, masking of f1 to zero is not performed.  It
% is performed by outside this function; that way the symmetry of the matrix can
% be more easily achieved.  The symmetry for the cases where f1 = 2 should be
% achieved separately, though.
Ds_cell = cell(1, v.count);
if s == Sign.p  % Ds == Df
	for w = v.elems
		f1 = 1;

		if bc(w) == BC.p
			fg = exp(-1i * gridnd.kBloch(w) * gridnd.L(w));
		else  % ghost point BC is BC.e if ge == GT.prim; BC.m if ge == GT.dual
			fg = 0;
		end

		Ds_cell{w} = create_Dw(w, N, f1, fg);
	end
else  % Ds == Db
	for w = v.elems
		if (ge == GT.prim && bc(w) == BC.m) || (ge == GT.dual && bc(w) == BC.e)
			f1 = 2;  % symmetry of operator for this case is not implemented yet
		else
			f1 = 1;
		end

		if bc(w) == BC.p
			fg = exp(-1i * gridnd.kBloch(w) * gridnd.L(w));
		else  % bc(w) == BC.e or BC.m
			fg = 0;  % f1 = 1 or 2 takes care of the ghost point
		end

		Ds_cell{w} = create_Dw(w, N, f1, fg);
		Ds_cell{w} = -Ds_cell{w}';  % conjugate transpose rather than transpose (hence nonsymmetry for kBloch ~= 0)
	end
end

dl = cell(1, v.count);
if isempty(dl_factor_cell)
	[dl{:}] = ndgrid(gridnd.dl{:,g});
else
	dl_cell = mult_vec(dl_factor_cell(:,g), gridnd.dl(:,g));
	[dl{:}] = ndgrid(dl_cell{:});
end

for w = v.elems
	Ds_cell{w} = create_spdiag(dl{w}.^-1) * Ds_cell{w};
end
