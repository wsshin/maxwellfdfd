function power = powerflux_patch(E_cell, H_cell, normal_axis, intercept, rect)

chkarg(istypesizeof(E_cell, 'Scalar3dcell', [1, Axis.count]), ...
	'"E_cell" should be length-%d row cell array with Scalar3d as elements.', Axis.count)
chkarg(istypesizeof(H_cell, 'Scalar3dcell', [1, Axis.count]), ...
	'"H_cell" should be length-%d row cell array with Scalar3d as elements.', Axis.count)
chkarg(istypesizeof(normal_axis, 'Axis'), '"normal_axis" should be instance of Axis.');
chkarg(istypesizeof(intercept, 'real'), '"intercept" should be real');

[h, v, n] = cycle(normal_axis);

Sn = poynting(n, E_cell{h}, E_cell{v}, H_cell{h}, H_cell{v}, n, intercept);

if nargin < 5  % no rect
	power = flux_patch(Sn);
else
	power = flux_patch(Sn, rect);
end

%%%
% Eh2d = slice_scalar3d(E_cell{h}, n, intercept);
% Ev2d = slice_scalar3d(E_cell{v}, n, intercept);
% Hh2d = slice_scalar3d(H_cell{h}, n, intercept);
% Hv2d = slice_scalar3d(H_cell{v}, n, intercept);
% 
% grid2d = Eh2d.grid2d;
% 
% hbound = rect(Dir.h,:);
% vbound = rect(Dir.v,:);
% chkarg(issorted(hbound) && issorted(vbound), 'each row of "rect" should be sorted in ascending order.');
% 
% chkarg(grid2d.contains(rect.'), '"rect" should be contained in grid.');
% 
% hall = grid2d.lall{Dir.h, GK.prim};
% vall = grid2d.lall{Dir.v, GK.prim};
% [Xh, Yv] = meshgrid(hall, vall);
% 
% in = find(hall > hbound(Sign.n), 1, 'first');
% ip = find(hall < hbound(Sign.p), 1, 'last');
% jn = find(vall > vbound(Sign.n), 1, 'first');
% jp = find(vall < vbound(Sign.p), 1, 'last');
% 
% hi = [hbound(Sign.n), hall(in:ip), hbound(Sign.p)];
% vi = [vbound(Sign.n), vall(jn:jp), vbound(Sign.p)];
% [XIh, YIv] = meshgrid(hi, vi);
% 
% Eh = Eh2d.array;
% Eh = permute(Eh, int([Dir.v, Dir.h]));
% Eh = interp2(Xh, Yv, Eh, XIh, YIv);
% Eh = ipermute(Eh, int([Dir.v, Dir.h]));
% 
% Ev = Ev2d.array;
% Ev = permute(Ev, int([Dir.v, Dir.h]));
% Ev = interp2(Xh, Yv, Ev, XIh, YIv);
% Ev = ipermute(Ev, int([Dir.v, Dir.h]));
% 
% Hh = Hh2d.array;
% Hh = permute(Hh, int([Dir.v, Dir.h]));
% Hh = interp2(Xh, Yv, Hh, XIh, YIv);
% Hh = ipermute(Hh, int([Dir.v, Dir.h]));
% 
% Hv = Hv2d.array;
% Hv = permute(Hv, int([Dir.v, Dir.h]));
% Hv = interp2(Xh, Yv, Hv, XIh, YIv);
% Hv = ipermute(Hv, int([Dir.v, Dir.h]));
% 
% Sn = real(Eh .* conj(Hv) - Ev .* conj(Hh)) / 2;
% 
% power = trapz(hi, Sn, int(Dir.h));  % integrate along horizontal direction
% power = trapz(vi, power);  % integrate along vertical direction
