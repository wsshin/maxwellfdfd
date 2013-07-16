function init_axes3d(axes_handle, grid3d, withinterp, withpml)
chkarg(~isempty(axes_handle) && ishandle(axes_handle), '"axes_handle" should be handle.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');
chkarg(istypesizeof(withinterp, 'logical'), '"withinterp" should be logical.');
chkarg(istypesizeof(withpml, 'logical'), '"withpml" should be logical.');

hold(axes_handle, 'on');

if withinterp
	lplot = grid3d.lplot(GT.prim, withinterp, withpml);
else
	lplot = grid3d.lpixelbound(GT.prim, withpml);
end
axis(axes_handle, [lplot{Axis.x}([1 end]), lplot{Axis.y}([1 end]), lplot{Axis.z}([1 end])]);
daspect(axes_handle, [1 1 1]);  % axis(ha, 'image');
box(axes_handle, 'on');
view(axes_handle, -38.5, 16)

hr = rotate3d;
set(hr, 'Enable', 'on');
setAllowAxesRotate(hr, axes_handle, true);

% set(axes_handle, 'FontSize', 18, 'FontWeight', 'Bold');
xlabel(axes_handle, char(Axis.x), 'Rotation', 0);
ylabel(axes_handle, char(Axis.y), 'Rotation', 0);
zlabel(axes_handle, char(Axis.z), 'Rotation', 0);
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
