function Ms_cell = create_Ms(s, ge, gridnd)
% Creates the forward or backward averaging matrix.  The output matrices perform
% simple, arithmetic averaging between two adjacent field components.

chkarg(istypesizeof(s, 'Sign'), '"s" should be instance of Sign.');  % Ds = Df (or Db) when s == Sign.p (or Sign.n)
chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for E-field
chkarg(istypesizeof(gridnd, 'Grid2d') || istypesizeof(gridnd, 'Grid3d'), '"gridnd" should be instance of Grid2d or Grid3d.');

v = Axis.x;
if istypesizeof(gridnd, 'Grid2d')
	v = Dir.h;
end

%% Get the shape.
N = gridnd.N;

%% Get the relevant derivative matrices.
bc = gridnd.bc;

% Basic setup of Mf and Mb.  Here, masking of f1 to zero is not performed.  It
% is performed by outside this function; that way the symmetry of the matrix can
% be more easily achieved, because both sides of each curl operator can be
% masked (see create_curl()).  The symmetry for the cases where f1 = 2 should be
% achieved separately, though.
isaddition = true;
Ms_cell = cell(1, v.count);
if s == Sign.p  % Ds == Df
	for w = v.elems
		f1 = 1;

		if bc(w) == BC.p
			fg = exp(-1i * gridnd.kBloch(w) * gridnd.L(w));
		else  % ghost point BC is BC.e if ge == GT.prim; BC.m if ge == GT.dual
			fg = 0;
		end

		Ms_cell{w} = create_Dw(w, N, f1, fg, isaddition);
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

		Ms_cell{w} = create_Dw(w, N, f1, fg, isaddition);
		Ms_cell{w} = Ms_cell{w}';  % conjugate transpose rather than transpose (hence nonsymmetry for kBloch ~= 0)
	end
end

for w = v.elems
	Ms_cell{w} = 0.5 .* Ms_cell{w};
end
