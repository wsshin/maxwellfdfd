%% powerflux_patch
% Calculate the power flux through a rectangular patch.

%%% Syntax
%  power = powerflux_patch(E_cell, H_cell, normal_axis, intercept)
%  power = powerflux_patch(E_cell, H_cell, normal_axis, intercept, rect)

%%% Parameters
% *Input*
%
% * |E_cell|: _E_-field in the format of |[Ex, Ey, Ez]|, where each |Ew| is an
% instance of <Scalar3d.html |Scalar3d|>.
% * |H_cell|: _H_-field in the format of |[Hx, Hy, Hz]|, where each |Hw| is an
% instance of <Scalar3d.html |Scalar3d|>.
% * |normal_axis|: axis normal to a rectangular patch.  It should be one of
% |Axis.x|, |Axis.y|, |Axis.z|.
% * |intercept|: location of the rectangular patch in the |normal_axis|
% direction.
% * |rect|: bounds of the rectangle in the plane.  For |normal_axis = Axis.y|,
% it is in the format of |[zmin zmax; xmin xmax]|.  If unassigned, the entire
% cross section normal to |normal_axis| is used.
%
% *Output*
% 
% * |power|: calculated power flux through the rectangular patch.

%%% Example
%   [E, H] = maxwell_run({ARGUMENTS});
%   power = powerflux_patch(E, H, Axis.z, 0, [0 200; 0 100]);

%%% See Also
% <powerflux_box.html |powerflux_box|>, <maxwell_run.html |maxwell_run|>

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
% hall = grid2d.lall{Dir.h, GT.prim};
% vall = grid2d.lall{Dir.v, GT.prim};
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
