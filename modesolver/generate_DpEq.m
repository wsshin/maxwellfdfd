function DpEq = generate_DpEq(p, bc, dp_prim)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_prim);
N = Nh*Nv;
DpEq = sparse(N,N);

% Every Eq component participates in two p-derivatives.
% In one curl loop it is substracted (-Eq), and in the other loop it is added (+Eq).

% -Eq
negEq = -ones(Nh,Nv);

% On p = 0 plane
%
% Do not handle bc(Pp, Sign.n) == BC.p here.  That case is handled by bc(Pp,
% Sign.p) == BC.p. Note that bc(Pp, Sign.n) == BC.p implies bc(Pp, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.Et0
    if p == Dir.h
        negEq(1,:) = 0;
    else 
        assert(p == Dir.v);
        negEq(:,1) = 0;
    end
end

% On q = 0 plane
if bc(Qq, Sign.n) == BC.En0
    if Qq == Dir.h
        negEq(:,1) = 0;
    else 
        assert(Qq == Dir.v);
        negEq(1,:) = 0;
    end
end


DpEq = spdiags(negEq(:), 0, DpEq);  % 0 means -Eq(h,v) is used to calculate Hr(h+-0,v).

% +Eq
posEq = ones(Nh,Nv);

% On p = Np plane
if bc(p, Sign.p) == BC.Et0 || bc(p, Sign.p) == BC.p  % Note that actually there are no other cases, since bc(Pp, Sign.p) cannot be BC.En0.
    if p == Dir.h
        posEq(1,:) = 0;
    else
        assert(p == Dir.v)
        posEq(:,1) = 0;
    end
end

if p == Dir.h
    DpEq = spdiags(posEq(:), 1, DpEq);  % 1 means +Eq(h,v) is used to calculate Hr(h-1,v).
else
    assert(p == Dir.v)
    DpEq = spdiags(posEq(:), Nh, DpEq);  % Nh means +Eq(h,v) is used to calculate Hr(h,v-1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpEq(j*Nh,(j-1)*Nh+1) = DpEq(j*Nh,(j-1)*Nh+1) + 1;
        end
    else
        for i = 1:Nh
            DpEq((Nv-1)*Nh+i,i) = DpEq((Nv-1)*Nh+i,i) + 1;
        end
    end
end

DpEq = spdiags(dp_prim(:), 0, N, N) \ DpEq;  % To calculate the curl, Ep shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
