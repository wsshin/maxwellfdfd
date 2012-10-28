function DpHr = generate_DpHr(p, bc, dp_dual)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_dual);
N = Nh*Nv;
DpHr = sparse(N,N);

% Every Hr component participates in two p-derivatives.
% In one curl loop it is substracted (-Hr), and in the other loop it is added (+Hr).

% -Hr
negHr = -ones(Nh,Nv);

% On p = 0 plane
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by bc(p,
% Sign.p) == BC.p. Note that bc(p, Sign.n) == BC.p implies bc(p, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.m
    if p == Dir.h
        negHr(1,:) = 0;
    else 
        assert(p == Dir.v);
        negHr(:,1) = 0;
    end
end

DpHr = spdiags(negHr(:), 0, DpHr);  % 0 means -Hr(x,y) is used to calculate Eq(x+-0,y).

% +Hr
posHr = ones(Nh,Nv);

% On p = Np plane
% Hr in the first cell in the p-axis is not used for +Hr.
if bc(p, Sign.p) == BC.m || bc(p, Sign.p) == BC.p  % Note that actually there are no other cases: bc(p, Sign.p) cannot be BC.e.
    if p == Dir.h
        posHr(1,:) = 0;
    else
        assert(p == Dir.v)
        posHr(:,1) = 0;
    end
end

if p == Dir.h
    DpHr = spdiags(posHr(:), 1, DpHr);  % 1 means +Hr(x,y) is used to calculate Eq(x-1,y).
else
    assert(p == Dir.v)
    DpHr = spdiags(posHr(:), Nh, DpHr);  % Nh means +Hr(x,y) is used to calculate Eq(x,y-1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpHr(j*Nh,(j-1)*Nh+1) = DpHr(j*Nh,(j-1)*Nh+1) + 1;
        end
    else
        for i = 1:Nh
            DpHr((Nv-1)*Nh+i,i) = DpHr((Nv-1)*Nh+i,i) + 1;
        end
    end
end

DpHr = spdiags(dp_dual(:), 0, N, N) \ DpHr;  % To calculate the curl, Hr shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
