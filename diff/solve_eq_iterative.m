function [E, H, relres, iter, resvec] = solve_eq_iterative(maxit, tol, F0type, eqtype, pml, omega, eps_cell, mu_cell, s_factor_cell, J_cell, M_cell, grid3d)

% [~, b, g_from_f, Op] = create_eq(eqtype, pml, omega, eps_cell, mu_cell, s_factor_cell, J_cell, M_cell, grid3d);
eq = MatrixEquation(eqtype, pml, omega, eps_cell, mu_cell, s_factor_cell, J_cell, M_cell, grid3d);
[Op, b, GfromF] = eq.matrixfree_op();  % Op has different ordering from A, so EH_from_fg() cannot be used
[A, b] = eq.matrix_op();
Op = A;

chkarg(isequal(F0type, 'zero') || isequal(F0type, 'rand'), '"F0type" should be either ''zero'' or ''rand''.');

if isequal(F0type, 'zero')
    x0 = zeros(size(b));
else  % F0type == 'rand'
    x0 = rand(size(b));
end

[f, ~, relres, iter, resvec] = bicg(Op, b, tol, maxit, [], [], x0);
save('resvec', 'resvec');

%%
if length(eqtype.f) == 1
	g = GfromF(f);
	
	if eqtype.f == FT.e
		e = f;
		h = g;
	else  % eqtype.f == FT.h
		e = g;
		h = f;
	end	
else  % use both E and H
	assert(isequal(eqtype.f, FT.elems));
	f = reshape(f, [], 2);
	e = f(:,1);
	h = f(:,2);
end

[E, H] = EH_from_eh(e, h, eqtype, grid3d);

% Test symmetry with respect to the plane bisecting the x-axis.
if false
	E = test_sym(Axis.x, A, b, E);
end


function E = test_sym(w, A, b, E)
Nw = size(E{Axis.x}, int(w));

Ex = E{Axis.x}; Ey = E{Axis.y}; Ez = E{Axis.z};
e = [Ex(:), Ey(:), Ez(:)];
e = e.';
e = e(:);
fprintf('norm(b - A*e)/norm(b) = %e\n', norm(b - A*e)/norm(b));

if mod(Nw, 2) == 0
	ind2 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind2{v}{w} = [1:Nw/2+1, Nw/2:-1:2];
		else
			ind2{v}{w} = [1:Nw/2, Nw/2:-1:1];
		end
	end
	
	ind3 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind3{v}{w} = [1, Nw:-1:Nw/2+1, Nw/2+2:Nw];
		else
			ind3{v}{w} = [Nw:-1:Nw/2+1, Nw/2+1:Nw];
		end
	end
else
	Nw_half = (Nw+1)/2;
	
	ind2 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind2{v}{w} = [1:Nw_half, Nw_half:-1:2];
		else
			ind2{v}{w} = [1:Nw_half, Nw_half-1:-1:1];
		end
	end
	
	ind3 = {{':', ':', ':'}, {':', ':', ':'}, {':', ':', ':'}};
	for v = Axis.elems
		if v == w
			ind3{v}{w} = [1, Nw:-1:Nw_half+1, Nw_half+1:Nw];
		else
			ind3{v}{w} = [Nw:-1:Nw_half, Nw_half+1:Nw];
		end
	end
end
Ex2 = Ex(ind2{Axis.x}{:});
Ey2 = Ey(ind2{Axis.y}{:});
Ez2 = Ez(ind2{Axis.z}{:});
E2 = {Ex2, Ey2, Ez2};

Ex3 = Ex(ind3{Axis.x}{:});
Ey3 = Ey(ind3{Axis.y}{:});
Ez3 = Ez(ind3{Axis.z}{:});
E3 = {Ex3, Ey3, Ez3};

Ex4 = (Ex2 + Ex3)/2;
Ey4 = (Ey2 + Ey3)/2;
Ez4 = (Ez2 + Ez3)/2;
E4 = {Ex4, Ey4, Ez4};

e2 = [Ex2(:), Ey2(:), Ez2(:)];
e2 = e2.';
e2 = e2(:);

e3 = [Ex3(:), Ey3(:), Ez3(:)];
e3 = e3.';
e3 = e3(:);

e4 = [Ex4(:), Ey4(:), Ez4(:)];
e4 = e4.';
e4 = e4(:);

fprintf('norm(b - A*e2)/norm(b) = %e\n', norm(b - A*e2)/norm(b));
fprintf('norm(b - A*e3)/norm(b) = %e\n', norm(b - A*e3)/norm(b));
fprintf('norm(b - A*e4)/norm(b) = %e\n', norm(b - A*e4)/norm(b));
fprintf('norm(b) = %e\n', norm(b));
fprintf('norm(e) = %e\n', norm(e));
fprintf('norm(A, 1) = %e\n', norm(A, 1));

E = E3;
