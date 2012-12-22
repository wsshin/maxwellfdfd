function ind = ind_for_loc(loc, axis, gk, grid3d)

ind = ismembc2(loc, grid3d.l{axis,gk});
if ind == 0
	[~, ind] = min(abs(grid3d.l{axis,gk} - loc));
	warning('FDS:srcAssign', ...
		['%s grid in %s-axis of "grid3d" does not have location %s; ', ...
		'closest grid vertex at %e will be taken instead.'], ...
		char(gk), char(axis), num2str(loc), grid3d.l{axis,gk}(ind));
end
