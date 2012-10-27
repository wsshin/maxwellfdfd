function lprim_part_cell = generate_lprim1d_part(domain, Lpml, interval_array, lprim0_array, ldual0_array)

% Check "domain".
chkarg(istypesizeof(domain, 'Interval'), '"domain" should be instance of Interval.');
for s = Sign.elems
	chkarg(domain.dl_boundary(s) == domain.dl_max, ...
		'"domain" should have the same dl_boundary as dl_max.');
end
dl_max = domain.dl_max;

% Check "Lpml".
chkarg(istypeof(Lpml, 'real'), 'element of "Lpml" should be real.');
chkarg(all(Lpml >= 0), 'element of "Lpml" should be positive.');
chkarg(isexpandable2row(Lpml, Sign.count), ...
	'"Lpml" should be scalar or length-%d vector.', Sign.count);
Lpml = expand2row(Lpml, Sign.count);
chkarg(domain.L >= Lpml(Sign.n) + Lpml(Sign.p), ...
	'size of "domain" should be at least the total depth of PML.');

% Check "intervals".
chkarg(istypesizeof(interval_array, 'Interval', [1 0]), ...
	'"interval_array" should be row vector of instances of Interval.');

% Check "lprim_array" and "ldual_array".
chkarg(istypesizeof(lprim0_array, 'real', [1 0]), '"lprim0_array" should be row vector with real elements.');
chkarg(istypesizeof(ldual0_array, 'real', [1 0]), '"ldual0_array" should be row vector with real elements.');

% Collect boundary points.
b_min = domain.bound(Sign.n);
b_max = domain.bound(Sign.p);
b_pml_min = b_min + Lpml(Sign.n);
b_pml_max = b_max - Lpml(Sign.p);
lprim0_array = [lprim0_array, b_min, b_pml_min, b_pml_max, b_max];  % possible duplicates in case of Lpml == 0 are removed later
	
% for i = 1:length(interval_array)
% 	for s = Sign.elems
% 		val = interval_array(i).bound(s);
% 		if val > b_min && val < b_max  % different from domain.contains(val)
% 			ldual0_array = [ldual0_array(1:end), val];
% 		end
% 	end
% end
bound = NaN(length(interval_array), Sign.count);
for i = 1:length(interval_array)
	bound(i,:) = interval_array(i).bound;
end
bound = reshape(bound, 1, []);
is_in = (bound > b_min) & (bound < b_max);
ldual0_array = [ldual0_array, bound(is_in)];


lprim0_array = unique(lprim0_array, 'sorted');  % duplicate elements are removed
assert(lprim0_array(1) == b_min && lprim0_array(end) == b_max);

ldual0_array = unique(ldual0_array, 'sorted');  % duplicate elements are removed
common = intersect(lprim0_array, ldual0_array);
assert(isempty(common), 'primary and dual grid share %s.', mat2str(common));

[l_array, ind] = sort([lprim0_array, ldual0_array]);
g_array = [GK.prim(ones(size(lprim0_array))), GK.dual(ones(size(ldual0_array)))];
g_array = g_array(ind);
n = length(l_array);

% For each boundary point, find the smallest dl suggested by intervals.
dl_array = NaN(1,n);
for j = 1:n
	val = l_array(j);
	
% 	% Handle intervals narrower than intervals' dl.
% 	if j == 1
% 		dl_dist = l_array(j+1) - val;
% 		g_nbr = g_array(j+1);
% 	elseif j == n
% 		dl_dist = val - l_array(j-1);
% 		g_nbr = g_array(j-1);
% 	else
% 		[dl_dist, ind] = min([val-l_array(j-1), l_array(j+1)-val]);
% 		g_nbr = g_array([j-1, j+1]);
% 		g_nbr = g_nbr(ind);
% 	end
% 	
% 	if g_array(j) ~= g_nbr  % consider cases where val and its neighbor are grid points with different GK
% 		dl_dist = 2/3 * dl_dist;
% 	end
% 	dl_min = min([dl_max, dl_dist]);

	dl_min = dl_max;
	
	for i = 1:length(interval_array)
		inter = interval_array(i);
		if inter.contains(val)
			if val == inter.bound(Sign.n)
				dl = inter.dl_boundary(Sign.n);
			elseif val == inter.bound(Sign.p)
				dl = inter.dl_boundary(Sign.p);
			else
				dl = inter.dl_max;
			end

			if dl < dl_min  % false if dl == NaN
				dl_min = dl;
			end
		end
	end
	
	dl_array(j) = dl_min;
end

% Create subgrids.
lprim_part_cell = {};
for j = 1:n
	val = l_array(j);
	dl = dl_array(j);
	if val == b_min
		lprim_part = [val, val+dl];
	elseif val == b_max
		lprim_part = [val-dl, val];
	elseif ismembc2(val, lprim0_array) ~= 0  % val is element of lprim0_array
		lprim_part = [val-dl, val, val+dl];
	else  % val is element of ldual0_array
		lprim_part = [val - dl/2, val + dl/2];
	end
	
	if isempty(lprim_part_cell)
		lprim_part_cell = {lprim_part};
	else
		lprim_part_cell_last = lprim_part_cell{end};
		
% 		if ~isempty(setdiff(lprim_part, lprim_part_cell_last))  % lprim_part is not contained in lprim_part_cell{end}
% 			lprim_last = lprim_part_cell_last(end);
% 			if lprim_part(1) > lprim_last
% 				lprim_part_cell = [lprim_part_cell(1:end), {lprim_part}];
% 			elseif lprim_part(1) == lprim_last
% 				lprim_part_cell{end} = unique([lprim_part_cell_last, lprim_part], 'sorted');
% 			else
% 				assert(lprim_part(1) < lprim_last);
% 				exception = MException('FDS:gridGen', 'subgrid generation failed: subgrids %s and %s overlap', ...
% 					mat2str(lprim_part_cell_last), mat2str(lprim_part));
% 				throw(exception);			
% 			end
% 		end

		common = intersect(lprim_part, lprim_part_cell_last);  % common is sorted
		if isempty(common)
			assert(lprim_part(1) > lprim_part_cell_last(end));
			lprim_part_cell = [lprim_part_cell(1:end), {lprim_part}];
		else
			n_overlap_curr = ismembc2(common(end), lprim_part);
			i_in_last = ismembc2(common(1), lprim_part_cell_last);
			n_overlap_last = length(lprim_part_cell_last) - i_in_last + 1;
			n_common = length(common);
			if n_overlap_curr == n_common && n_overlap_last == n_common
				lprim_part_cell{end} = [lprim_part_cell_last(1:i_in_last-1), lprim_part];
			else
				exception = MException('FDS:gridGen', 'subgrid generation failed: subgrids %s and %s overlap but cannot be combined', ...
					mat2str(lprim_part_cell_last), mat2str(lprim_part));
				throw(exception);
			end
		end			
	end
end

% if all(lprim_part_cell{1} == lprim_part_cell{end})  % initial and final subgrids are the same
% 	lprim_part_cell = lprim_part_cell(1);  % == {lprim_part_cell{1}}
% end
