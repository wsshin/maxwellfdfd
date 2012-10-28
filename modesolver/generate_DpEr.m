function DpEr = generate_DpEr(p, bc, dp_prim)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_prim);
N = Nh*Nv;
DpEr = sparse(N,N);

% Every Er component participates in two x-derivatives.
% In one curl loop it is substracted (-Er), and in the other loop it is added (+Er).

% +Er
posEr = ones(Nh,Nv);  

% On p = 0 plane
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by
% bc(p,Sign.p) == BC.p. Note that bc(p,Sign.n) == BC.p implies bc(p,Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.e
    if p == Dir.h
        posEr(1,:) = 2;
    else
        assert(p == Dir.v);
        posEr(:,1) = 2;
    end
end

if bc(p, Sign.n) == BC.m
    if p == Dir.h
        posEr(1,:) = 0;
    else
        assert(p == Dir.v);
        posEr(:,1) = 0;
    end
end

DpEr = spdiags(posEr(:), 0, DpEr);  % 0 means +Er(x,y) is used to calculate Hq(x,y).

% -Er
negEr = -ones(Nh,Nv);  

% On p = Np plane
% Er in the last cell in the p-axis is not used for -Er.
if p == Dir.h
    negEr(Nh,:) = 0;
    DpEr = spdiags(negEr(:), -1, DpEr);  % -1 means -Er(x,y) is used to calculate Hq(x+1,y).
else
    assert(p == Dir.v);
    negEr(:,Nv) = 0;
    DpEr = spdiags(negEr(:), -Nh, DpEr);  % -Nh means -Er(x,y) is used to calculate Hq(x,y+1).
end 

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 0:Nv-1
            DpEr(j*Nh+1,(j+1)*Nh) = DpEr(j*Nh+1,(j+1)*Nh) - 1;
        end
    else
        assert(p == Dir.v);
        for i = 1:Nh
            DpEr(i,(Nv-1)*Nh+i) = DpEr(i,(Nv-1)*Nh+i) - 1;
        end
    end
end

DpEr = spdiags(dp_prim(:), 0, N, N) \ DpEr;  % To calculate the curl, Er shuold be divided by the edge lengths centered at primary grid vertices.
