function power = powerflux_box(E_cell, H_cell, box)

chkarg(istypesizeof(E_cell, 'Scalar3dcell', [1, Axis.count]), ...
	'"E_cell" should be length-%d row cell array with Scalar3d as elements.', Axis.count)
chkarg(istypesizeof(H_cell, 'Scalar3dcell', [1, Axis.count]), ...
	'"H_cell" should be length-%d row cell array with Scalar3d as elements.', Axis.count)
chkarg(istypesizeof(box, 'real', [Axis.count, Sign.count]), ...
	'"box" should be %d-by-%d array with real elements.', Axis.count, Sign.count);
chkarg(issorted(box(Axis.x,:)) && issorted(box(Axis.y,:)) && issorted(box(Axis.z,:)), ...
	'each row of "box" should be sorted in ascennding order.');

p = NaN(Axis.count, Sign.count);
for n = Axis.elems
	[h, v] = cycle(n);
	for s = Sign.elems
		intercept = box(n,s);
		rect = box([h v], :);
		p(n,s) = powerflux_patch(E_cell, H_cell, n, intercept, rect);
		if s == Sign.n
			p(n,s) = -p(n,s);  % power flux going out of box
		end
	end
end

power = sum(p(:));

