function format_colorbar(axes_handle, scalar)
chkarg(istypesizeof(scalar, 'Scalar2d') || istypesizeof(scalar, 'Scalar3d'), ...
	'"scalar" shuold be instancef of either Scalar2d or Scalar3d.');

% axes_handle = colorbar;

% Make sure that 0 is ticked when 0 is the minimum of the colorbar.
crange = caxis(gca);  % caxis is the property of the axes handle for a figure
ticks = get(axes_handle, 'YTick');
if crange(1) == 0 && ticks(1) > 0
	crange(1) = -eps * crange(2);
	caxis(gca, crange);
	ticks = [0 ticks];
end
set(axes_handle, 'YTick', ticks);

% Prevent 'x10^n' from showing on top of the colorbar.
ticklabel = textscan(num2str(ticks),'%s');
ticklabel = ticklabel{1};
set(axes_handle, 'YTickLabel', ticklabel);

% if scalar.unitvalue ~= 1e0
% 	str = [' (', char(hex2dec('00D7')), num2str(scalar.unitvalue, '%.2e'), ')'];
% 	title(axes_handle, str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');
% end
str = [' (', char(hex2dec('00D7'))];
if scalar.unitvalue < 1e4 && scalar.unitvalue > 1e-3
	str = [str, num2str(scalar.unitvalue), ')'];
else
	str = [str, num2str(scalar.unitvalue, '%.2e'), ')'];
end
title(axes_handle, str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');
