function [E, H] = EH_from_eh(e, h, eqtype, grid3d)

N = grid3d.N;
e = reshape(e, Axis.count, prod(N));
Ex = e(int(Axis.x), :); Ex = reshape(Ex, N);
Ey = e(int(Axis.y), :); Ey = reshape(Ey, N); 
Ez = e(int(Axis.z), :); Ez = reshape(Ez, N);
E = {Ex, Ey, Ez};

h = reshape(h, Axis.count, prod(N));
Hx = h(int(Axis.x), :); Hx = reshape(Hx, N);
Hy = h(int(Axis.y), :); Hy = reshape(Hy, N); 
Hz = h(int(Axis.z), :); Hz = reshape(Hz, N);
H = {Hx, Hy, Hz};