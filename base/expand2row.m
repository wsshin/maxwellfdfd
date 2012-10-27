function row = expand2row(obj, n)
if isscalar(obj)
	row = repmat(obj, 1, n);
else
	chkarg(isvector(obj) && length(obj)==n, '"obj" should be length-"n" vector.');
	if iscolumn(obj)
		obj = obj.';
	end
	row = obj;
end
