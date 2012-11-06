function collapsed_array = collapse_complex(ri_array)

% Note that if "complex_array" is a column vector, then its real and imaginary
% parts are the row vectors of "expanded_array", because the fastest-varying
% index (i.e., index in the column direction) is assigned for alternating the
% real and imaginary parts.
chkarg(istypeof(ri_array, 'real'), '"ri_array" shuold be array with real elements.');
nD = ndims(ri_array);
ind_r = cell(1, nD);
ind_i = cell(1, nD);
for d = 2:nD
	ind_r{d} = ':';
	ind_i{d} = ':';
end
ind_r{1} = 1;
ind_i{1} = 2;

collapsed_array = ri_array(ind_r{:}) + 1i * ri_array(ind_i{:});

% Permute array indices to, for example, (axis, x, y, z, ri); note that ri is
% always 1, so the indices are actually (axis, x, y, z).
collapsed_array = permute(collapsed_array, [2:nD, 1]);
