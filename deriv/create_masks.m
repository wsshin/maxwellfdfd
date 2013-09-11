function [ind_Mp, ind_Md] = create_masks(ge, gridnd)
% Generating the indices to mask is better than generating a diagonal mask
% matrix, because multiplying the mask matrix cannot make NaN or Inf elements
% zeros.

chkarg(istypesizeof(ge, 'GT'), '"ge" should be instance of GT.');  % ge: grid type for the E-field
chkarg(istypesizeof(gridnd, 'Grid2d') || istypesizeof(gridnd, 'Grid3d'), '"gridnd" should be instance of Grid2d or Grid3d.');

v = Axis.x;
if istypesizeof(gridnd, 'Grid2d')
	v = Dir.h;
end

%% Get the shape.
N = gridnd.N;
bc = gridnd.bc;
ind = cell(1, v.count);

% Mask matrices
mask_p = cell(1, v.count);
for w = v.elems  % mask w-component of field on primary grid
	mask = false(N);
	us = setdiff(v.elems, w);

	for u = us
		if (ge == GT.prim && bc(u) == BC.e) || (ge == GT.dual && bc(u) == BC.m)
			[ind{:}] = deal(':');  % ind = {':', ':', ':'} for grid3d
			ind{u} = 1;
			mask(ind{:}) = true;
		end
	end
	mask_p{w} = mask(:);
end

ind_Mp = cell2mat(mask_p);
ind_Mp = ind_Mp(:);

mask_d = cell(1, v.count);
for w = v.elems  % mask w-component of field on dual grid
	mask = false(N);
	
	if (ge == GT.prim && bc(w) == BC.e) || (ge == GT.dual && bc(w) == BC.m)
		[ind{:}] = deal(':');  % ind = {':', ':', ':'} for grid3d
		ind{w} = 1;
		mask(ind{:}) = true;
	end
	mask_d{w} = mask(:);
end

ind_Md = cell2mat(mask_d);
ind_Md = ind_Md(:);

