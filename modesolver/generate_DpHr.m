function DpHr = generate_DpHr(p, bc, dp_dual)

chkarg(istypesizeof(p, 'Dir'), '"p" should be instance of Dir.');

[Nh, Nv] = size(dp_dual);
N = Nh*Nv;
DpHr = sparse(N,N);

% Every Hr component participates in two x-derivatives.
% In one curl loop it is substracted (-Hr), and in the other loop it is added (+Hr).

% +Hr
posHr = ones(Nh,Nv);  

% On p = 0 plane
%
% Do not handle bc(Pp, Sign.n) == BC.p here.  That case is handled by bc(Pp,Pos)
% == BC.p. Note that bc(Pp, Sign.n) == BC.p implies bc(Pp,Pos) == BC.p.
if bc(p, Sign.n) == BC.En0
    if p == Dir.h
        posHr(1,:) = 2;
    else
        assert(p == Dir.v);
        posHr(:,1) = 2;
    end
end

DpHr = spdiags(posHr(:), 0, DpHr);  % 0 means +Hr(x,y) is used to calculate Er(x,y).

% -Hr
negHr = -ones(Nh,Nv);  

% On p = Np plane
if p == Dir.h
    negHr(Nh,:) = 0;
    DpHr = spdiags(negHr(:), -1, DpHr);  % -1 means -Hr(x,y) is used to calculate Er(x+1,y).
else
    assert(p == Dir.v);
    negHr(:,Nv) = 0;
    DpHr = spdiags(negHr(:), -Nh, DpHr);  % -Nx means -Hr(x,y) is used to calculate Ex(x,y+1).
end 

if bc(p, Sign.p) == BC.p
    if p == Dir.h
        for j = 0:Nv-1
            DpHr(j*Nh+1,(j+1)*Nh) = DpHr(j*Nh+1,(j+1)*Nh) - 1;
        end
    else
        assert(p == Dir.v);
        for i = 1:Nh
            DpHr(i,(Nv-1)*Nh+i) = DpHr(i,(Nv-1)*Nh+i) - 1;
        end
    end
end

DpHr = spdiags(dp_dual(:), 0, N, N) \ DpHr;  % To calculate the curl, Hr shuold be divided by the edge lengths centered at primary grid vertices.
