function expanded_array = expand_complex(complex_array)

% Note that if "complex_array" is a column vector, then its real and imaginary
% parts are the row vectors of "expanded_array", because the fastest-varying
% index (i.e., index in the column direction) is assigned for alternating the
% real and imaginary parts.
chkarg(istypeof(complex_array, 'complex'), '"complex_array" should be array with complex elements.');
nD = ndims(complex_array);
expanded_array = cat(nD+1, real(complex_array), imag(complex_array));
expanded_array = permute(expanded_array, [nD+1, 1:nD]);
