%% complete_lprim1d
% Complete a 1D primary grid by filling the gaps between parts of the grid.

%%% Syntax
%  lprim_array = complete_lprim1d(lprim_part_cell)

%%% Description
% |complete_lprim1d(lprim_part_cell)| attempts to generate a full 1D primary
% grid from a partially generated grid.  
%
% |lprim_part_cell| is a cell, whose element is either a number |dl| (grid cell
% size) or an array |[l_i, ..., l_f]| (_subgrid_, which is a complete grid
% between |l_i| and |l_f|).  A real number |dl|, which comes between two
% neighboring subgrids, describes the target grid cell size between two
% subgrids.
%
% An example of |lprim_part_cell| is
%
%  {[0 0.5], 1, [10 11 12], 1.5, [20.5, 21], 2, [30, 31]}
%
% where |[0.0 0.5]|, |[10 11 12]|, |[20.5, 21]|, |[30, 31]| are subgrids; |1| is
% the target |dl| between |0.5| and |10|; |1.5| is the target |dl| between |12|
% and |20.5|; |2| is the target |dl| between |21| and |30|.  Then,
% |complete_lprim1d(lprim_part_cell)| generates a grid between |0| and |31|.  An
% error is generated when the function fails to generate a primary grid.

%%% Example
%   % Complete a primary grid.
%   lprim_part_cell = {[0 0.5], 1, [10 11 12], 1.5, [20.5, 21], 2, [30, 31]};
%   lprim_array = complete_lprim1d(lprim_part_cell);
%
%   % Test the quality of the completed grid.
%   deltas = diff(grid);
%   ratios = deltas(1:end-1)./deltas(2:end);
%   is_smooth = all(ratios > 1/1.3) && all(ratios < 1.3);
%   fprintf('Does the completed primary grid have smoothly varying separations?  ');
%   if is_smooth
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end

function lprim_array = complete_lprim1d(lprim_part_cell)

% Check inputs.
chkarg(istypesizeof(lprim_part_cell, 'realcell', [1 0], [1 0]), '"lprim_part_cell" should be row cell array with real row vectors as elements.');
numelems = length(lprim_part_cell);
chkarg(mod(numelems, 2) == 1, '"lprim_part_cell" should have odd number of elements.');

ds = {};  % cell array whose element is {subgrid, dl_target}, except for last element being {subgrid}
for i = 1:(numelems+1)/2
	if i == (numelems+1)/2
		curr = lprim_part_cell(2*i-1);
	else
		curr = lprim_part_cell([2*i-1, 2*i]);
		chkarg(istypesizeof(curr{2}, 'real'), 'element #%d of "lprim_part_cell" should be real.', 2*i);
	end
	chkarg(length(curr{1})>=2 && issorted(curr{1}), ...
		'element #%d of "lprim_part_cell" should be length-2 or longer row vector with real elements in ascending order.', 2*i-1);

	if i >= 2
		chkarg(prev{1}(end) < curr{1}(1), 'subgrids %s and %s in "lprim_part_cell" should be sorted in ascending order.', mat2str(prev{1}), mat2str(curr{1}));
	end
	ds = {ds{1:end}, curr};
	prev = curr;
end

if numelems == 1
	lprim_array = ds{1}{1};
	return
end

% Fill gaps between neighboring subgrids.
rt = 1.9;  % target ratio of geometric sequence; rt = 1.2 to be more strict
rmax = 2.0;  % maximum ratio of geometric sequence; rmax = 1.3 to be more strict
numgrids = length(ds);
assert(numgrids >= 2);
curr = ds{1};  %  curr == {subgrid, dl_target} or curr == {subgrid}
lprim_array = curr{1};
for i = 2:numgrids
	dl_n = lprim_array(end) - lprim_array(end-1);  % == curr(end) - curr(end-1); dl of current subgrid
	next = ds{i};  %  next == {subgrid, dl_target} or next == {subgrid}
	dl_p = next{1}(2) - next{1}(1);  % dl of next subgrid
	gap = [lprim_array(end), next{1}(1)];  % == [curr(end), next{1}(1)]; gap between neighboring subgrids
	dl_target = curr{2};
	try
		filler = fill_targeted_geometric(dl_n, gap, dl_target, dl_p, rt, rmax);
	catch err
		exception = MException('Maxwell:gridGen', ['grid generation failed between subgrids ', ...
			mat2str(curr{1}), ' and ', mat2str(next{1}), ' with target dl = %s: %s'], ...
			num2str(curr{2}), err.message);
		throw(exception);
	end
	
	% Below, note that the provided subgrids are preserved in spite of roundoff
	% errors in "filler".  This is important for assiging materials and sources
	% at intended locations.  For example, failure to maintain the provided
	% subgrids can result in warnings in Shape.generate_kernel().
	lprim_array = [lprim_array(1:end), filler(2:end-1), next{1}];
	curr = next;
end

dlprim_array = diff(lprim_array);
ind = find_stiff_ddl(dlprim_array, rmax);
if ~isempty(ind)
	exception = MException('Maxwell:gridGen', ['grid generation failed: grid vertex locations ', ...
		'[..., %e, %e, %e, ...] are separated by [..., %e, %e, ...] that do not vary smoothly.'], ...
		lprim_array(ind(1)), lprim_array(ind(1)+1), lprim_array(ind(1)+2), dlprim_array(ind(1)), dlprim_array(ind(1)+1));
	throw(exception);
end


function ind = find_stiff_ddl(dl_array, rt)
chkarg(istypesizeof(dl_array, 'real', [1 0]) && length(dl_array) >= 2, ...
	'"dl_array" should be length-2 or longer row vector with real elements in ascending order.');
dl1 = dl_array(1:end-1);
dl2 = dl_array(2:end);

if rt >= 1
	ind = [find(dl1./dl2 > rt), find(dl1./dl2 < 1/rt)];
else
	ind = [find(dl1./dl2 > 1/rt), find(dl1./dl2 , rt)];
end


function truth = isequal_approx(a, b)
d = abs(a-b);
truth = d < min([a b]) * 1e-8;

function truth = is_smooth(dl_array, rt)
truth = isempty(find_stiff_ddl(dl_array, rt));


function filler = fill_constant(dl_min, dl_max, gap, rt, rmax)
chkarg(dl_min <= dl_max || isequal_approx(dl_min, dl_max), ...
	'"dl_min = %s" should not be greater than "dl_max = %s".', num2str(dl_min), num2str(dl_max));
chkarg(is_smooth([dl_min, dl_max], rt), '"dl_min" and "dl_max" should be similar.');

L = gap(2) - gap(1);
chkarg(L > 0, 'second element of "gap" should be greater than the first element.');

numcells = [floor(L/dl_max), ceil(L/dl_max)];  % candidate numbers of grid cells
dl = L ./ numcells;  % length(dl) == 2

assert(is_smooth([dl_min, dl(1)], rmax) || is_smooth([dl_min, dl(2)], rmax), ...
	'dl = %e is too small or dl = %e is too large for gap size = %e.', dl_min, dl_max, L);

% Unless the coarser grid fails to generate smoothly varying grid, gives priority to it.
if ~is_smooth([dl_min, dl(1)], rmax)
	numcells = numcells(2);  % finer grid
else
	numcells = numcells(1);  % coarser grid
end
% 	filler = gap(1):dl:gap(end);  % this may not include gap(end) due to roundoff error
filler = linspace(gap(1), gap(2), numcells+1);


function filler = fill_targeted_geometric_sym(dl_sym, gap, dl_t, rt, rmax)
% Generate a filler subgrid whose grid vertex separations increases
% geometrically from both ends.  The target dl should be equal to or greater
% than dl_n before the gap and dl_p after the gap, but dl_n == dl_p.

% gap: length-2 row vector whose elements are the locations of the beginning and
% ending grid vertices of the gap
% dl_sym: dl of the subgrids before and after the gap, i.e., dl_n == dl_p
% dl_t: target dl in the gap
% rt: target ratio of geometric sequence
% rmax: maximum ratio of geometric sequence

chkarg(dl_t >= dl_sym || isequal_approx(dl_t, dl_sym), ...
	'"dl_t = %s" should be equal to or greater than "dl_sym = %s".', num2str(dl_t), num2str(dl_sym));

L = gap(2) - gap(1);
chkarg(L > 0, 'second element of "gap" should be greater than the first element.');

if is_smooth([dl_sym, dl_t], rt)
	filler = fill_constant(dl_sym, dl_t, gap, rt, rmax);
else
	dl_max = dl_t;
	dl_min = dl_sym;
	
	% Guess graded dl's.
	n = ceil(log(dl_max/dl_min)/log(rt));  % smallest n satisfying (dl_max/dl_min)^(1/n) <= rt
	r = (dl_max/dl_min)^(1/n);  % ratio of geometric sequence
	dl_array = dl_min * (r.^(1:n));
	
	L_graded = sum(dl_array);  % dl_min * (r^1 + ... + r^n)
	if 2*L_graded > L
		% Try best to geometrically increase dl close to dl_t.
		n = NaN(1, 2);
		n(1) = fzero(@(n) dl_min*rt * (rt^(n-1) + 2*(rt^(n-1) - 1)/(rt-1)) - L, 1);  % ceil(n(1)) is the smallest integer satisfying dl_min * (rt^1 + ... + rt^n + ... + rt^1) >= L; note that rt^n is added once
		n(2) = fzero(@(n) dl_min*rt * 2*(rt^n - 1)/(rt-1) - L, 1);  % ceil(n(2)) is the smallest integer satisfying dl_min * (rt^1 + ... + rt^n + rt^n ... + rt^1) >= L; note that rt^n is added twice
		[~, i] = min(abs(ceil(n) - n));  % between n(1) and n(2), choose the one closer to their smallest integers bounding from above
		n = ceil(n(i));
		assert(n >= 1 , 'dl = %e is too small or dl = %e is too large for gap size = %e.', dl_min, dl_max, L);
		if i == 1
			r = fzero(@(s) dl_min*s * (s^(n-1) + 2*(s^(n-1) - 1)/(s-1)) - L, rt);  % dl_min * (r^1 + ... + r^n + ... + r^1) == L
			assert(r <= rt);
			assert(r >= 1/rt, 'dl = %e is too small or dl = %e is too large for gap size = %e.', dl_min, dl_max, L);
			dl_array = dl_min * (r.^[1:n, n-1:-1:1]);
		else
			assert(i==2);
			r = fzero(@(s) dl_min*s * (s^n - 1)/(s-1) - L/2, rt);  % dl_min * (r^1 + ... + r^n + r^n + ... + r^1) == L
			assert(r <= rt);
			assert(r >= 1/rt, 'dl = %e is too small or dl = %e is too large for gap size = %e.', dl_min, dl_max, L);
			dl_array = dl_min * (r.^[1:n, n:-1:1]);
		end			
		dl_filler = dl_array;
	else
		% Slightly under-fill the gap with dl_max and the above generated graded
		% dl's.
		n_dl_max = floor((L - 2*L_graded)/dl_max);  % use floor() to under-fill
		dl_max_array = dl_max(ones(1, n_dl_max));  % [dl_max, ..., dl_max]
		L_dl_max = sum(dl_max_array);

		% Slightly over-fill the gap with dl_min, graded dl's, and the above
		% generated L_dl_max.
		assert(2*L_graded + L_dl_max <= L);
		n_dl_min = ceil((L - 2*L_graded - L_dl_max)/2/dl_min);  % use ceil() to over-fill
		dl_min_array = dl_min(ones(1, n_dl_min));  % [dl_min, ..., dl_min]
		L_dl_min = sum(dl_min_array);

		% Update the graded dl's.
		L_graded = (L - L_dl_max - 2*L_dl_min)/2;
		r = fzero(@(s) dl_min*s * (s^n - 1) / (s-1) - L_graded, r);  % dl_min * (r^1 + ... + r^n) == L_graded
		dl_array = dl_min * (r.^(1:n));
		dl_filler = [dl_min_array, dl_array, dl_max_array, fliplr(dl_array), dl_min_array];		
	end
		
	% Construct filler.
	filler = cumsum([gap(1), dl_filler]);  % last element of filler is gap(2)
end


function filler = fill_targeted_geometric(dl_n, gap, dl_t, dl_p, rt, rmax)
% Generate a filler subgrid whose grid vertex separations increases
% geometrically from both ends.

% gap: length-2 row vector whose elements are the locations of the beginning and
% ending grid vertices of the gap
% dl_n, dl_p: dl of the subgrids before and after the gap
% dl_t: target dl in the gap
% rt: target ratio of geometric sequence
% rmax: maximum ratio of geometric sequence

chkarg((dl_t >= dl_n || isequal_approx(dl_t, dl_n)) && (dl_t >= dl_p || isequal_approx(dl_t, dl_p)), ...
	'"dl_t = %s" should be equal to or greater than "dl_n = %s" and "dl_p = %s".', ...
	num2str(dl_t), num2str(dl_n), num2str(dl_p));

L = gap(2) - gap(1);
chkarg(L > 0, 'second element of "gap" should be greater than the first element.');

if dl_n == dl_p
	filler = fill_targeted_geometric_sym(dl_n, gap, dl_t, rt, rmax);
else  % dl_n < dl_p or dl_n > dl_p
	dl_max = max(dl_n, dl_p);  
	dl_min = min(dl_n, dl_p);
	
	% Guess graded dl's.
	n = ceil(log(dl_max/dl_min)/log(rt));  % smallest n satisfying (dl_max/dl_min)^(1/n) <= rt
	r = (dl_max/dl_min)^(1/n);  % ratio of geometric sequence
	dl_array = dl_min * (r.^(1:n));
	
	% Slightly under-fill the gap with dl_max and the above generated graded
	% dl's.
	L_graded = sum(dl_array);  % dl_min * (r^1 + ... + r^n)
	assert(L_graded <= L, 'dl = %e is too small or dl = %e is too large for gap size = %e.', dl_min, dl_max, L);
	if dl_n < dl_p
		filler_n = cumsum([gap(1), dl_array]);
		gap_sym = [filler_n(end), gap(2)];
		filler_sym = fill_targeted_geometric_sym(dl_p, gap_sym, dl_t, rt, rmax);
		filler = [filler_n(1:end-1), filler_sym];
	else
		assert(dl_n > dl_p)
		filler_p = cumsum([gap(2), -dl_array]);
		filler_p = fliplr(filler_p);
		gap_sym = [gap(1), filler_p(1)];
		filler_sym = fill_targeted_geometric_sym(dl_n, gap_sym, dl_t, rt, rmax);
		filler = [filler_sym, filler_p(2:end)];
	end	
end
