function [truth, normal_axis] = is2dlike(N)

chkarg(istypesizeof(N, 'int', [1, Axis.count]), ...
	'"N" should be length-%d row vector with integer elements.');

[~, imin] = min(N);
n = Axis.elems(imin);
[h, v] = cycle(n);

truth =  N(n)/N(h) < 0.1 && N(n)/N(v) < 0.1;
normal_axis = n;
