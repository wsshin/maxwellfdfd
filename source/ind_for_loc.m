function ind = ind_for_loc(loc, axis, gt, grid3d)

ind = ismembc2(loc, grid3d.l{axis,gt});
if ind == 0
	[~, ind] = min(abs(grid3d.l{axis,gt} - loc));
	warning('Maxwell:srcAssign', ...
		['%s grid in %s-axis of "grid3d" does not have location %s; ', ...
		'closest grid vertex at %e will be taken instead.'], ...
		char(gt), char(axis), num2str(loc), grid3d.l{axis,gt}(ind));
end
