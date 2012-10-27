function array = cell2array(cellarray, nD)

% The reason for not taking nD as ndims(cellarray{1}) is that an array
% intended to be 3D can look 2D because the size of the 3rd dimension is 1.
chkarg(istypesizeof(nD, 'int') && nD > 0, '"nD" should be positive integer.');
chkarg(istypesizeof(cellarray, 'complexcell', [1 0], zeros(1, nD)), ...
	'"cellarray" should be a row cell array whose each element is %dD array with complex elements.', nD);

ncell = length(cellarray);
if ncell == 0
	array = [];
else
	assert(ncell >= 1);
	dims = size(cellarray{1});
	chkarg(istypesizeof(cellarray, 'complexcell', [1 0], dims), 'elements of "cellarray" should be have the same size.');
	array = cat(nD+1, cellarray{:});
	array = permute(array, [nD+1, 1:nD]);
end
