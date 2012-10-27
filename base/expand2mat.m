function mat = expand2mat(obj, nrow, ncol)

if isscalar(obj)
	mat = repmat(obj, nrow, ncol);
elseif isvector(obj)
	chkarg(length(obj)==nrow, 'length of "obj" array should be "nrow".');
	if isrow(obj)
		obj = obj.';
	end
	mat = repmat(obj, 1, ncol);
else
	chkarg(ismatrix(obj) && all(size(obj)==[nrow, ncol]), ...
		'size of "obj" should be "nrow"-by-"ncol".')
	mat = obj;
end
