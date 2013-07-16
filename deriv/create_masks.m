function [Mp_cell, Md_cell] = create_masks(ge, grid3d)

chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for the E-field
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

%% Get the shape.
N = grid3d.N;
bc = grid3d.bc;

my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));

% Mask matrices
mask_p = cell(1, Axis.count);
for w = Axis.elems  % mask w-component of Fp
	mask = ones(N);
	[u, v] = cycle(w);

	if (ge == GT.prim && bc(u) == BC.e) || (ge == GT.dual && bc(u) == BC.m)
		ind = {':', ':', ':'};
		ind{u} = 1;
		mask(ind{:}) = 0;
	end
	if (ge == GT.prim && bc(v) == BC.e) || (ge == GT.dual && bc(v) == BC.m)
		ind = {':', ':', ':'};
		ind{v} = 1;
		mask(ind{:}) = 0;
	end
	mask_p{w} = mask;
end
Mp_cell = my_diag([mask_p{Axis.x}(:); mask_p{Axis.y}(:); mask_p{Axis.z}(:)]);

mask_d = cell(1, Axis.count);
for w = Axis.elems  % mask w-component of Fd
	mask = ones(N);
	if (ge == GT.prim && bc(w) == BC.e) || (ge == GT.dual && bc(w) == BC.m)
		ind = {':', ':', ':'};
		ind{w} = 1;
		mask(ind{:}) = 0;
	end
	mask_d{w} = mask;
end
Md_cell = my_diag([mask_d{Axis.x}(:); mask_d{Axis.y}(:); mask_d{Axis.z}(:)]);

