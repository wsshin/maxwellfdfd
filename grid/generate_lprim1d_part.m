function lprim_part_cell = generate_lprim1d_part(domain, Lpml, interval_array, lprim0_array, ldual0_array)

% Check "domain".
chkarg(istypesizeof(domain, 'Interval'), '"domain" should be instance of Interval.');
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
lprim_inter = [];
for i = 1:length(interval_array)
	lprim_inter = [lprim_inter(1:end), interval_array(i).lprim];
end
is_in = (lprim_inter > b_min) & (lprim_inter < b_max);
lprim0_array = [lprim0_array, lprim_inter(is_in)];  % boundaries of shapes are aligned with primary gird

lprim0_array = unique(lprim0_array);  % sorted and duplicate elements are removed
isequal_approx = @(a, b) abs(a-b) < dl_max * 1e-8;
dlprim0_array = diff(lprim0_array);
ind_unique = ~isequal_approx(dlprim0_array, 0);
lprim0_array = [lprim0_array(ind_unique), lprim0_array(end)];
if isequal_approx(lprim0_array(end), lprim0_array(end-1))
	lprim0_array = lprim0_array(1:end-1);
end

assert(lprim0_array(1) == b_min && lprim0_array(end) == b_max);

ldual0_array = unique(ldual0_array);  % sorted and duplicate elements are removed
common = intersect(lprim0_array, ldual0_array);
if ~isempty(common)
	exception = MException('Maxwell:gridGen', 'primary and dual grid share %s.', mat2str(common));
	throw(exception);
end

% [l_array, ind] = sort([lprim0_array, ldual0_array]);
% g_array = [GT.prim(ones(size(lprim0_array))), GT.dual(ones(size(ldual0_array)))];
% g_array = g_array(ind);
% n = length(l_array);s


% Generates primary grid points around dual grid points.
dl_prim0_array = diff(lprim0_array);
dl_dual0_array = diff(ldual0_array);
dl_dual0_array = min([dl_dual0_array, Inf; Inf, dl_dual0_array]);  % min along columns

n_dual0 = length(ldual0_array);
lprim_by_ldual0 = [];
for j = 1:n_dual0
	val = ldual0_array(j);
	ind = find(lprim0_array < val, 1, 'last');
	dl_min = min([dl_max, dl_dual0_array(j), dl_prim0_array(ind)]);
	
	for i = 1:length(interval_array)
		inter = interval_array(i);
		if inter.contains(val)
			dl = inter.dl_max;
			if dl < dl_min
				dl_min = dl;
			end
		end
	end
	
	lprim_new = val - dl_min/2;
	if lprim_new >= b_min && lprim_new <= b_max
		lprim_by_ldual0 = [lprim_by_ldual0(1:end), lprim_new];
	end
	
	lprim_new = val + dl_min/2;
	if lprim_new >= b_min && lprim_new <= b_max
		lprim_by_ldual0 = [lprim_by_ldual0(1:end), lprim_new];
	end
end

lprim0_array = [lprim0_array, lprim_by_ldual0];
lprim0_array = unique(lprim0_array);  % sorted and duplicate elements are removed

% For each interval between primary grid points, find the smallest dl suggested by intervals.
lprim0_mid_array = (lprim0_array(1:end-1) + lprim0_array(2:end)) / 2;
dl_prim0_array = diff(lprim0_array);

n_prim0 = length(lprim0_array);
dl_mid_array = NaN(1, n_prim0-1);
for j = 1:n_prim0-1
	val = lprim0_mid_array(j);
	dl_min = min([dl_max, dl_prim0_array(j)]);
	
	for i = 1:length(interval_array)
		inter = interval_array(i);
		if inter.contains(val)
			dl = inter.dl_max;
			if dl < dl_min
				dl_min = dl;
			end
		end
	end
	
	dl_mid_array(j) = dl_min;
end

dl_boundary_array = min([Inf, dl_mid_array; dl_mid_array, Inf]);  % min along columns

% Create subgrids.  (This is the part that needs to be improved.)
val = lprim0_array(1);
dl = dl_boundary_array(1);
prev = [val, val+dl];
lprim_part_cell = {prev};
for j = 2:n_prim0
	val = lprim0_array(j);
	dl = dl_boundary_array(j);
	if val == b_max
		curr = [val-dl, val];
	else  % val ~= b_min or b_max
		curr = [val-dl, val, val+dl];  %  same cell size on both sides of primary node; important for eps = 2/(1/eps1 + 1/eps2)
	end
	
	if isequal_approx(curr(1), prev(end-1)) && isequal_approx(curr(2), prev(end))  % curr(1) == prev(end-1) && curr(2) == prev(end)
		curr = [prev, curr(3:end)];
		lprim_part_cell = [lprim_part_cell(1:end-1), {curr}];
	elseif isequal_approx(curr(1), prev(end))  % curr(1) == prev(end)
		curr = [prev, curr(2:end)];
		lprim_part_cell = [lprim_part_cell(1:end-1), {curr}];
	else
		lprim_part_cell = [lprim_part_cell(1:end), {dl_mid_array(j-1), curr}];
	end
	
	prev = curr;
end
