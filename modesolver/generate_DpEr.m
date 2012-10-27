function DpEr = generate_DpEr(p, bc, dp_prim)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_prim);
N = Nh*Nv;
DpEr = sparse(N,N);

% Every Er component participates in two p-derivatives.
% In one curl loop it is substracted (-Er), and in the other loop it is added (+Er).

% -Er
negEr = -ones(Nh,Nv);

% On p = 0 plane
%
% Do not handle bc(Pp, Sign.n) == BC.p here.  That case is handled by bc(Pp,
% Sign.p) == BC.p. Note that bc(Pp, Sign.n) == BC.p implies bc(Pp, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.Et0
    if p == Dir.h
        negEr(1,:) = 0;
    else 
        assert(p == Dir.v);
        negEr(:,1) = 0;
    end
end

DpEr = spdiags(negEr(:), 0, DpEr);  % 0 means -Er(x,y) is used to calculate Hr(x+-0,y).

% +Er
posEr = ones(Nh,Nv);

% On p = Np plane
if bc(p, Sign.p) == BC.Et0 || bc(p, Sign.p) == BC.p  % Note that actually there are no other cases, since bc(Pp, Sign.p) cannot be BC.En0.
    if p == Dir.h
        posEr(1,:) = 0;
    else
        assert(p == Dir.v)
        posEr(:,1) = 0;
    end
end

if p == Dir.h
    DpEr = spdiags(posEr(:), 1, DpEr);  % 1 means +Er(x,y) is used to calculate Hr(x-1,y).
else
    assert(p == Dir.v)
    DpEr = spdiags(posEr(:), Nh, DpEr);  % Nx means +Er(x,y) is used to calculate Hr(x,y-1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpEr(j*Nh,(j-1)*Nh+1) = DpEr(j*Nh,(j-1)*Nh+1) + 1;
        end
    else
        for i = 1:Nh
            DpEr((Nv-1)*Nh+i,i) = DpEr((Nv-1)*Nh+i,i) + 1;
        end
    end
end

DpEr = spdiags(dp_prim(:), 0, N, N) \ DpEr;  % To calculate the curl, Ep shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
