function DpHp = generate_DpHp(p, bc, dp_prim)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_prim);
N = Nh*Nv;
DpHp = sparse(N,N);

% Every Hp component participates in two p-derivatives.
% In one derivative it is substracted (-Hp), and in the other derivative it is added (+Hp).

% +Hp
posHp = ones(Nh,Nv);

% BC.m boundary condition on p = 0 plane
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by bc(p,
% Sign.p) == BC.p. Note that bc(p, Sign.n) == BC.p implies bc(p, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.m  % Hp is symmetric around p = 0, so subtraction becomes zero
    if p == Dir.h
        posHp(1,:) = 0;
    else 
        assert(p == Dir.v);
        posHp(:,1) = 0;
    end
end

if bc(p, Sign.n) == BC.e  % Hp is antisymmetric around p = 0
    if p == Dir.h
        posHp(1,:) = 2;
    else 
        assert(p == Dir.v);
        posHp(:,1) = 2;
    end
end

DpHp = spdiags(posHp(:), 0, DpHp);  % 0 means +Hp(h,v) is used to calculate rho(h,v).

% -Hp
negHp = -ones(Nh,Nv);

% BC.m and periodic boundary condition on p = Np plane
% Hp in the last cell in the p-axis is not used for -Hp.
if bc(p, Sign.p) == BC.m || bc(p, Sign.p) == BC.p  % Note that actually there are no other cases: bc(p, Sign.p) cannot be BC.e.
    if p == Dir.h
        negHp(end,:) = 0;
    else
        assert(p == Dir.v)
        negHp(:,end) = 0;
    end
end

if p == Dir.h
    DpHp = spdiags(negHp(:), -1, DpHp);  % -1 means -Hp(h,v) is used to calculate rho(h+1,v).
else
    assert(p == Dir.v)
    DpHp = spdiags(negHp(:), -Nh, DpHp);  % -Nh means -Hp(h,v) is used to calculate rho(h,v+1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpHp((j-1)*Nh+1,j*Nh) = DpHp((j-1)*Nh+1,j*Nh) - 1;
        end
    else
        for i = 1:Nh
            DpHp(i,(Nv-1)*Nh+i) = DpHp(i,(Nv-1)*Nh+i) - 1;
        end
    end
end

DpHp = spdiags(dp_prim(:), 0, N, N) \ DpHp;  % To calculate the curl, Hp shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
