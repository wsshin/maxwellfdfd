function phi = flux_patch(scalar2d, rect)

chkarg(istypesizeof(scalar2d, 'Scalar2d'), '"scalar2d" should be instance of Scalar2d.');

grid2d = scalar2d.grid2d;
if nargin < 2  % no rect
	rect = grid2d.bound_plot(false);  % false: do not include PML
end
chkarg(istypesizeof(rect, 'real', [2 2]), '"rect" should be %d-by-%d array with real elements.');

chkarg(issorted(rect(Dir.h,:)) && issorted(rect(Dir.v,:)), 'each row of "rect" should be sorted in ascending order.');
chkarg(grid2d.contains(rect(Dir.h,:), rect(Dir.v,:)), '"rect" should be contained in grid.');

gt = scalar2d.gt_array;
% l = grid2d.lall(Dir.elems + Dir.count*subsindex(gt));
% la = grid2d.lall(Dir.elems + Dir.count*subsindex(alter(gt)));

la = grid2d.lvoxelbound(gt, true);  % true: include PML
array = scalar2d.data_original();

ind_array = NaN(Dir.count, Sign.count);
dl = cell(1, Dir.count);
for d = Dir.elems
	ind_dl = NaN(1, Sign.count);
	for s = Sign.elems
		bound = rect(d, s);
		i = ismembc2(bound, la{d});
		if i == 0  % bound is not in l{d}
			warning('Maxwell:interp', ['patch boundary %s = %s is not aligned with %s grid;', ...
				'flux will be slightly inaccurate.'], ...
				char(grid2d.axis(d)), num2str(bound), char(alter(gt(d))));
		end
			
		if s == Sign.n
			ind_dl(s) = find(la{d} > bound, 1, 'first');
			ind_array(d,s) = ind_dl(s) - 1;
		else  % s == Sign.p
			ind_dl(s) = find(la{d} < bound, 1, 'last');
			ind_array(d,s) = ind_dl(s);
		end
	end
	dl{d} = [rect(d,Sign.n), la{d}(ind_dl(Sign.n):ind_dl(Sign.p)), rect(d,Sign.p)];
	dl{d} = diff(dl{d});
end

% Calculate the Rieman sum.
array = array(...
	ind_array(Dir.h,Sign.n):ind_array(Dir.h,Sign.p), ...
	ind_array(Dir.v,Sign.n):ind_array(Dir.v,Sign.p));

phi = dl{Dir.h} * array * dl{Dir.v}.';
