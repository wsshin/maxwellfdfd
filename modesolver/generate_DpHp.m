function DpHp = generate_DpHp(p, bc, dp_prim)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_prim);
N = Nh*Nv;
DpHp = sparse(N,N);

% Every Hp component participates in two p-derivatives.
% In one derivative it is substracted (-Hp), and in the other derivative it is added (+Hp).

% -Hp
negHp = -ones(Nh,Nv);

% BC.Et0 boundary condition on p = 0 plane
%
% Do not handle bc(Pp, Sign.n) == BC.p here.  That case is handled by bc(Pp,
% Sign.p) == BC.p. Note that bc(Pp, Sign.n) == BC.p implies bc(Pp, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.Et0
    if p == Dir.h
        negHp(1,:) = 0;
    else 
        assert(p == Dir.v);
        negHp(:,1) = 0;
    end
end

DpHp = spdiags(negHp(:), 0, DpHp);  % 0 means -Hp(x,y) is used to calculate Hr(x+-0,y).

% +Hp
posHp = ones(Nh,Nv);

% BC.Et0 and periodic boundary condition on p = Np plane
if bc(p, Sign.p) == BC.Et0 || bc(p, Sign.p) == BC.p  % Note that actually there are no other cases, since bc(Pp, Sign.p) cannot be BC.En0.
    if p == Dir.h
        posHp(1,:) = 0;
    else
        assert(p == Dir.v)
        posHp(:,1) = 0;
    end
end

if p == Dir.h
    DpHp = spdiags(posHp(:), 1, DpHp);  % 1 means +Hp(x,y) is used to calculate Hr(x-1,y).
else
    assert(p == Dir.v)
    DpHp = spdiags(posHp(:), Nh, DpHp);  % Nx means +Hp(x,y) is used to calculate Hr(x,y-1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpHp(j*Nh,(j-1)*Nh+1) = DpHp(j*Nh,(j-1)*Nh+1) + 1;
        end
    else
        for i = 1:Nh
            DpHp((Nv-1)*Nh+i,i) = DpHp((Nv-1)*Nh+i,i) + 1;
        end
    end
end

DpHp = spdiags(dp_prim(:), 0, N, N) \ DpHp;  % To calculate the curl, Ep shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
