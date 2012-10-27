function DpHq = generate_DpHq(p, bc, dp_dual)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

q = p.alter();

[Nh, Nv] = size(dp_dual);
N = Nh*Nv;
DpHq = sparse(N,N);

% Every Hq component participates in two p-derivatives. In one curl loop it is
% substracted (-Hq), and in the other loop it is added (+Hq).

% +Hq
posHq = ones(Nh,Nv);  

% BC.Et0 boundary condition on q = 0 plane: Ep = 0 for q = 0 => Hq = 0.
%
% We don't need to handle the periodic or BC.En0 boundary condition in the
% q-direction.
if bc(q, Sign.n) == BC.Et0
    if q == Dir.h  % p == Dir.v
        posHq(1,:) = 0;
    else  % p == Dir.h
        assert(q == Dir.v);
        posHq(:,1) = 0;
    end
end

% BC.En0 boundary condition on p = 0 plane: can assume Hq = 0 for p = 0.
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by bc(p,
% Sign.p) == BC.p. Note that bc(p, Sign.n) == BC.p implies bc(p, Sign.p)
% == BC.p.
if bc(p, Sign.n) == BC.En0
    if p == Dir.h
        posHq(1,:) = 2;
    else
        assert(p == Dir.v);
        posHq(:,1) = 2;
    end
end

% BC.Et0 boundary condition on p = 0 plane: can assume the same Hq behind the p =
% 0 plane.
if bc(p, Sign.n) == BC.Et0
    if p == Dir.h
        posHq(1,:) = 0;
    else
        assert(p == Dir.v);
        posHq(:,1) = 0;
    end
end

DpHq = spdiags(posHq(:), 0, DpHq);  % 0 means +Hq(h,v) is used to calculate Er(h,v).

% -Hq
negHq = -ones(Nh,Nv);  

% BC.Et0 boundary condition on q = 0 plane
%
% We don't need to handle the periodic or BC.En0 boundary condition in the
% q-direction.
if bc(q, Sign.n) == BC.Et0
    if q == Dir.h  % p == Dir.v
        negHq(1,:) = 0;
    else  % p == Dir.h
        assert(q == Dir.v);
        negHq(:,1) = 0;
    end
end

% On p = Lp plane
if p == Dir.h
    negHq(Nh,:) = 0;
    DpHq = spdiags(negHq(:), -1, DpHq);  % -1 means -Hq(h,v) is used to calculate Er(h+1,v).
else
    assert(p == Dir.v);
    negHq(:,Nv) = 0;
    DpHq = spdiags(negHq(:), -Nh, DpHq);  % -Nh means -Hq(h,v) is used to calculate Er(h,v+1).
end 

% Periodic boundary condition on p = Lp plane
if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 0:Nv-1
            DpHq(j*Nh+1,(j+1)*Nh) = DpHq(j*Nh+1,(j+1)*Nh) - 1;
        end
    else
        assert(p == Dir.v);
        for i = 1:Nh
            DpHq(i,(Nv-1)*Nh+i) = DpHq(i,(Nv-1)*Nh+i) - 1;
        end
    end
end

DpHq = spdiags(dp_dual(:), 0, N, N) \ DpHq;  % To calculate the curl, Hq shuold be divided by the edge lengths centered at primary grid vertices.
