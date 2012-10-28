function DpEq = generate_DpEq(p, bc, dp_prim)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');
q = alter(p);

[Nh, Nv] = size(dp_prim);
N = Nh*Nv;
DpEq = sparse(N,N);

% Every Eq component participates in two p-derivatives. In one curl loop it is
% substracted (-Eq), and in the other loop it is added (+Eq).

% +Eq
posEq = ones(Nh,Nv);  

% BC.m boundary condition on q = 0 plane: Ep = 0 for q = 0 => Eq = 0.
%
% We don't need to handle the periodic or BC.e boundary condition in the
% q-direction.
if bc(q, Sign.n) == BC.m
    if q == Dir.h  % p == Dir.v
        posEq(1,:) = 0;
    else  % p == Dir.h
        assert(q == Dir.v);
        posEq(:,1) = 0;
    end
end

% BC.e boundary condition on p = 0 plane: can assume Eq = 0 for p = 0.
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by bc(p,
% Sign.p) == BC.p. Note that bc(p, Sign.n) == BC.p implies bc(p, Sign.p)
% == BC.p.
if bc(p, Sign.n) == BC.e
    if p == Dir.h
        posEq(1,:) = 2;
    else
        assert(p == Dir.v);
        posEq(:,1) = 2;
    end
end

% BC.m boundary condition on p = 0 plane: can assume the same Eq behind the p =
% 0 plane.
if bc(p, Sign.n) == BC.m
    if p == Dir.h
        posEq(1,:) = 0;
    else
        assert(p == Dir.v);
        posEq(:,1) = 0;
    end
end

DpEq = spdiags(posEq(:), 0, DpEq);  % 0 means +Eq(h,v) is used to calculate Hr(h,v).

% -Eq
negEq = -ones(Nh,Nv);  

% BC.m boundary condition on q = 0 plane
%
% We don't need to handle the periodic or BC.e boundary condition in the
% q-direction.
if bc(q, Sign.n) == BC.m
    if q == Dir.h  % p == Dir.v
        negEq(1,:) = 0;
    else  % p == Dir.h
        assert(q == Dir.v);
        negEq(:,1) = 0;
    end
end

% On p = Lp plane
% Eq in the last cell in the p-axis is not used for -Eq.
if p == Dir.h
    negEq(Nh,:) = 0;
    DpEq = spdiags(negEq(:), -1, DpEq);  % -1 means -Eq(h,v) is used to calculate Hr(h+1,v).
else
    assert(p == Dir.v);
    negEq(:,Nv) = 0;
    DpEq = spdiags(negEq(:), -Nh, DpEq);  % -Nh means -Eq(h,v) is used to calculate Hr(h,v+1).
end 

% Periodic boundary condition on p = Lp plane
if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 0:Nv-1
            DpEq(j*Nh+1,(j+1)*Nh) = DpEq(j*Nh+1,(j+1)*Nh) - 1;
        end
    else
        assert(p == Dir.v);
        for i = 1:Nh
            DpEq(i,(Nv-1)*Nh+i) = DpEq(i,(Nv-1)*Nh+i) - 1;
        end
    end
end

DpEq = spdiags(dp_prim(:), 0, N, N) \ DpEq;  % To calculate the curl, Eq shuold be divided by the edge lengths centered at primary grid vertices.
