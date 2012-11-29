function array = flip_array(array)

nD = ndims(array);
for i = 1:nD
	array = flipdim(array, i);
end
