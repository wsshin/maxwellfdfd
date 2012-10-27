function truth = isexpandable2row(obj, n)

truth = isscalar(obj) || (isvector(obj) && length(obj)==n);
