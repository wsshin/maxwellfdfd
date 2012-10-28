function DpHq = generate_DpHq(p, bc, dp_dual)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');
q = alter(p);

[Nh, Nv] = size(dp_dual);
N = Nh*Nv;
DpHq = sparse(N,N);

% Every Hq component participates in two p-derivatives.
% In one curl loop it is substracted (-Hq), and in the other loop it is added (+Hq).

% -Hq
negHq = -ones(Nh,Nv);

% On p = 0 plane
%
% Do not handle bc(p, Sign.n) == BC.p here.  That case is handled by bc(p,
% Sign.p) == BC.p. Note that bc(p, Sign.n) == BC.p implies bc(p, Sign.p) ==
% BC.p.
if bc(p, Sign.n) == BC.m
    if p == Dir.h
        negHq(1,:) = 0;
    else 
        assert(p == Dir.v);
        negHq(:,1) = 0;
    end
end

% % On q = 0 plane
% if bc(q, Sign.n) == BC.e
%     if q == Dir.h
%         negHq(:,1) = 0;
%     else 
%         assert(q == Dir.v);
%         negHq(1,:) = 0;
%     end
% end


DpHq = spdiags(negHq(:), 0, DpHq);  % 0 means -Hq(h,v) is used to calculate Er(h+-0,v).

% +Hq
posHq = ones(Nh,Nv);

% On p = Np plane
% Hq in the first cell in the p-axis is not used for +Hq.
if bc(p, Sign.p) == BC.m || bc(p, Sign.p) == BC.p  % Note that actually there are no other case: bc(p,Sign.p) cannot be BC.e.
    if p == Dir.h
        posHq(1,:) = 0;
    else
        assert(p == Dir.v)
        posHq(:,1) = 0;
    end
end

if p == Dir.h
    DpHq = spdiags(posHq(:), 1, DpHq);  % 1 means +Hq(h,v) is used to calculate Er(h-1,v).
else
    assert(p == Dir.v)
    DpHq = spdiags(posHq(:), Nh, DpHq);  % Nh means +Hq(h,v) is used to calculate Er(h,v-1).
end
    

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 1:Nv
            DpHq(j*Nh,(j-1)*Nh+1) = DpHq(j*Nh,(j-1)*Nh+1) + 1;
        end
    else
        for i = 1:Nh
            DpHq((Nv-1)*Nh+i,i) = DpHq((Nv-1)*Nh+i,i) + 1;
        end
    end
end

DpHq = spdiags(dp_dual(:), 0, N, N) \ DpHq;  % To calculate the curl, Hq shuold be divided by the primary edge lengths, which are centered at dual grid vertices.
