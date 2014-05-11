function init_axes2d(axes_handle, grid2d, withinterp, withpml, isswapped)
chkarg(~isempty(axes_handle) && ishandle(axes_handle), '"axes_handle" should be handle.');
chkarg(istypesizeof(grid2d, 'Grid2d'), '"grid2d" should be instance of Grid2d.');
chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');
chkarg(istypesizeof(isswapped, 'logical'), '"isswapped" should be logical.');

hold(axes_handle, 'on');

if withinterp
	lplot = grid2d.lplot(GT.prim, withinterp, withpml);
else
	lplot = grid2d.lvoxelbound(GT.prim, withpml);
end

if ~isswapped
	h = Dir.h; v = Dir.v;
else
	h = Dir.v; v = Dir.h;
end

axis(axes_handle, [lplot{h}([1 end]), lplot{v}([1 end])]);
daspect(axes_handle, [1 1 1]);  % axis(axes_handle, 'image');
box(axes_handle, 'on');

hr = rotate3d;
setAllowAxesRotate(hr, axes_handle, false);

% set(axes_handle, 'FontSize', 18, 'FontWeight', 'Bold');
str = [char(grid2d.axis(h)), '  (', char(hex2dec('00D7'))];
if grid2d.unitvalue < 1e5 && grid2d.unitvalue > 1e-3
	str = [str, num2str(grid2d.unitvalue), ')'];
else
	str = [str, num2str(grid2d.unitvalue, '%.2e'), ')'];
end

xlabel(axes_handle, str, 'Rotation', 0, 'FontSize', 15); 
% s = [char(grid2d.axis(v)), '  (', char(hex2dec('00D7')), num2str(grid2d.unitvalue, '%e'), ')'];
str = char(grid2d.axis(v));
ylabel(axes_handle, str, 'Rotation', 0, 'FontSize', 15); 

set(axes_handle, 'TickDir', 'out');


% L0 = scalar3d.gi.L0;
% dlu = scalar3d.gi.display_length_unit;
% dlu_scale = scalar3d.gi.dlu_scale;
% xlabel(strcat(AxisName(Xx), ' (', char(hex2dec('00D7')), ...  % 00D7 is a unicode character for the multiplication sign.
% 	num2str(L0*dlu_scale), dlu, ')'));
% ylabel(strcat(AxisName(Yy), ' (', char(hex2dec('00D7')), ...  % 00D7 is a unicode character for the multiplication sign.
% 	num2str(L0*dlu_scale), dlu, ')'));
% zlabel(strcat(AxisName(Zz), ' (', char(hex2dec('00D7')), ...  % 00D7 is a unicode character for the multiplication sign.
% 	num2str(L0*dlu_scale), dlu, ')'));
% title(plot_title);
