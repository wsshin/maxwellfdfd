%% istypesizeof
% Check the type and size of an object.

%%% Syntax
%  truth = istypesizeof(obj, typename)
%  truth = istypesizeof(obj, typename, dims)
%  truth = istypesizeof(obj, typename, dims, dimselem)

%%% Description
% |istypesizeof(obj, typename)| returns |true| if |obj| is a single scalar of
% the type described by a string |typename|.  It is equivalent to
% |istypesizeof(obj, typename, [1 1])|.  Exmaples of |typename| are |'int'|,
% |'real'|, |'complex'|, |'logical'| (for 'true' and 'false'), and the names of classes.
%
% |istypesizeof(obj, typename, dims)| returns |true| if |obj| is an array of the
% type described by |typename| and |size(obj) == dims|.
%
% |dims| must be a row vector of integers.  Some elements of |dims| can be |0|.
% In such cases, |istypesizeof(obj, typename, dims)| returns |true| no matter
% how many elements are in the dimensions associated with |0| in |dims|.  For
% example, |istypesizeof(obj, typename, [1 0])| returns |true| if |obj| is a row
% vector of any length, and |istypesizeof(obj, typename, [0 1])| returns |true|
% if |obj| is a column vector of any length.
%
% To test if |obj| is a cell array whose elements are scalars or arrays of a
% specific type, a _cell type_ can be specified as |typename|.  A _cell type_ is
% a conventional type name followed by 'cell', such as |'intcell'|,
% |'realcell'|, |'complexcell'|.
%
% |istypesizeof(obj, typename, dims, dimselem)| has an additional parameter
% |dimselem| that is effective only when |typename| is a cell type.  |dimselem|
% specifies the size of each element of the cell. For example,
% |istypesizeof(obj, 'realcell', [2 3], [1 0])| tests if |obj| is a 2-by-3 cell
% array whose each element is a row vector with real elements.
%
% For |obj == []|, the function returns |true| for any |typename| that is a
% non-cell type and any |dims|. For |obj == {}|, the function returns |true| for
% any |typename| that is a cell type and any |dims|.

%%% Example
%   % Test if obj is a scalar.
%   obj = 3+4i; type = 'complex';
%   fprintf('Test 1: is the object a %s scalar?  ', type);
%   if istypesizeof(obj, type)
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end
%
%   % Test if obj is an array of a specific size.
%   obj = rand(17, 17); type = 'real';
%   [m, n] = size(obj);
%   fprintf('Test 2: is the object a %s array of size %d-by-%d?  ', type, m, n);
%   if istypesizeof(obj, type, [n n]);
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end
%
%   % Test if obj is a cell array whose elements are scalars or arrays of Axis.
%   obj = {Axis.x, [Axis.y, Axis.z]; Axis.y, [Axis.z Axis.y; Axis.x, Axis.y]}; type = 'Axiscell';
%   [m, n] = size(obj);
%   fprintf('Test 3: is the object a cell of size %d-by-%d whose elements are scalars or arrays of %s?  ', m, n, type);
%   if istypesizeof(obj, type, [m n]);
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end
%
%   % Test if obj is a cell array whose elements are real arrays of a specific size.
%   n = 3;
%   elem1 = rand(n, n); elem2 = rand(n, n); elem3 = rand(n, n); elem4 = rand(n, n);
%   obj = {elem1, elem2; elem3, elem4}; type = 'realcell';
%   [m, n] = size(obj);
%   [me, ne] = size(elem1);
%   fprintf('Test 4: is the object a cell of size %d-by-%d whose elements are %s arrays of size %d-by-%d?  ', m, n, type, me, ne);
%   if istypesizeof(obj, type, [m n], [me ne]);
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end

%%% See Also
% <istypeof.html istypeof>

function truth = istypesizeof(obj, typename, dims, dimselem)

chkarg(ischar(typename), '"typename" should be string.');
iscelltype = false;
if ~isempty(strfind(typename, 'cell'))
	typename = strrep(typename, 'cell', '');
	iscelltype = true;
end

if nargin < 3  % no dims
	dims = [1 1];
end
chkarg(isrow(dims) && istypeof(dims, 'int'), ...
	'"dims" should be row vector with integer elements.');

size_obj = size(obj);
if length(size_obj) < length(dims)
	% An m-dimensional tensor can be seen as an n-dimensional tensor as well for n >= m.
	temp = size_obj;
	size_obj = ones(1, length(dims));
	size_obj(1:length(temp)) = temp;
end

if isempty(obj) && ~all(dims)  % if some element of "dims" is zero, allow empty "obj"
	truth = true;
elseif length(size_obj) ~= length(dims)
	truth = false;
else
	assert(length(size_obj) == length(dims));
	ind_zero = find(dims == 0);
	dims(ind_zero) = size_obj(ind_zero);

	if ~all(size_obj == dims)
		truth = false;
	elseif iscelltype
		truth = iscell(obj);
		obj = obj(:);		
		if nargin < 4  % no dimselem
			for i = 1:length(obj)
				if ~truth
					break
				end
				truth = truth && istypeof(obj{i}, typename);
			end
		else
			for i = 1:length(obj)
				if ~truth
					break
				end
				truth = truth && istypesizeof(obj{i}, typename, dimselem);
			end			
		end
	else
		truth = istypeof(obj, typename);
	end
end
