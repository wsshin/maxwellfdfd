function truth = isexpandable2mat(obj, nrow, ncol)

truth = isscalar(obj) || (isvector(obj) && length(obj)==nrow) || (ismatrix(obj) && all(size(obj)==[nrow, ncol]));
