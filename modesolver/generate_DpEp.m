function DpEp = generate_DpEp(p, bc, dp_prim)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_prim);
N = Nh*Nv;
DpEp = sparse(N,N);

% Every Ep component participates in two p-derivatives.
% In one derivative it is substracted (-Ep), and in the other derivative it is added (+Ep).

% -Ep
negEp = -ones(Nh,Nv);

% BC.m boundary condition on p = 0 plane
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by bc(p,
% Sign.p) == BC.p. Note that bc(p, Sign.n) == BC.p implies bc(p, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.m
    if p == Dir.h
        negEp(1,:) = 0;
    else 
        assert(p == Dir.v);
        negEp(:,1) = 0;
    end
end

DpEp = spdiags(negEp(:), 0, DpEp);  % 0 means -Ep(x,y) is used to calculate Er(x+-0,y).

% +Ep
posEp = ones(Nh,Nv);

% BC.m and periodic boundary condition on p = Np plane
% Ep in the first cell in the p-axis is not used for +Ep.
if bc(p, Sign.p) == BC.m || bc(p, Sign.p) == BC.p  % Note that actually there are no other cases: bc(p, Sign.p) cannot be BC.e.
    if p == Dir.h
        posEp(1,:) = 0;
    else
        assert(p == Dir.v)
        posEp(:,1) = 0;
    end
end

if p == Dir.h
    DpEp = spdiags(posEp(:), 1, DpEp);  % 1 means +Ep(x,y) is used to calculate Er(x-1,y).
else
    assert(p == Dir.v)
    DpEp = spdiags(posEp(:), Nh, DpEp);  % Nh means +Ep(x,y) is used to calculate Er(x,y-1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpEp(j*Nh,(j-1)*Nh+1) = DpEp(j*Nh,(j-1)*Nh+1) + 1;
        end
    else
        for i = 1:Nh
            DpEp((Nv-1)*Nh+i,i) = DpEp((Nv-1)*Nh+i,i) + 1;
        end
    end
end

DpEp = spdiags(dp_prim(:), 0, N, N) \ DpEp;  % To calculate the curl, Ep shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
