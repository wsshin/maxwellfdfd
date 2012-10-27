function DpEp = generate_DpEp(p, bc, dp_dual)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_dual);
N = Nh*Nv;
DpEp = sparse(N,N);

% Every Ep component participates in two p-derivatives.
% In one derivative it is substracted (-Ep), and in the other derivative it is added (+Ep).

% +Ep
posEp = ones(Nh,Nv);

% BC.Et0 boundary condition on p = 0 plane
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by bc(p,
% Sign.p) == BC.p. Note that bc(p, Sign.n) == BC.p implies bc(p, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.Et0  % Ep is symmetric around p = 0
    if p == Dir.h
        posEp(1,:) = 0;
    else 
        assert(p == Dir.v);
        posEp(:,1) = 0;
    end
end

if bc(p, Sign.n) == BC.En0  % Ep is antisymmetric around p = 0
    if p == Dir.h
        posEp(1,:) = 2;
    else 
        assert(p == Dir.v);
        posEp(:,1) = 2;
    end
end

DpEp = spdiags(posEp(:), 0, DpEp);  % 0 means +Ep(h,v) is used to calculate rho(h,v).

% -Ep
negEp = -ones(Nh,Nv);

% BC.Et0 and periodic boundary condition on p = Np plane
if bc(p, Sign.p) == BC.Et0 || bc(p, Sign.p) == BC.p  % Note that actually there are no other cases, since bc(p, Sign.p) cannot be BC.En0.
    if p == Dir.h
        negEp(end,:) = 0;
    else
        assert(p == Dir.v)
        negEp(:,end) = 0;
    end
end

if p == Dir.h
    DpEp = spdiags(negEp(:), -1, DpEp);  % -1 means -Ep(h,v) is used to calculate rho(h+1,v).
else
    assert(p == Dir.v)
    DpEp = spdiags(negEp(:), -Nh, DpEp);  % -Nh means -Ep(h,v) is used to calculate rho(h,v+1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpEp((j-1)*Nh+1,j*Nh) = DpEp((j-1)*Nh+1,j*Nh) - 1;
        end
    else
        for i = 1:Nh
            DpEp(i,(Nv-1)*Nh+i) = DpEp(i,(Nv-1)*Nh+i) - 1;
        end
    end
end

DpEp = spdiags(dp_dual(:), 0, N, N) \ DpEp;  % To calculate the curl, Ep shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
