%% istypeof
% Check the type of an object.

%%% Syntax
%  truth = istypeof(obj, typename)

%%% Description
% |istypeof(obj, typename)| returns |true| if |obj| is a scalar or array of the
% type described by a string |typename|. Exmaples of |typename| are |'int'|,
% |'real'|, |'complex'|, |'logical'| (for 'true' and 'false'), and the names of classes.  
%
% To test if |obj| is a cell array whose elements are scalars or arrays of a
% specific type, a _cell type_ can be specified as |typename|.  A _cell_type_ is
% a conventional type name followed by 'cell', such as |'intcell'|,
% |'realcell'|, |'complexcell'|.
%
% For |obj == []|, the function returns |true| for any |typename| that is a
% non-cell type. For |obj == {}|, the function returns |true| for any |typename|
% that is a cell type.

%%% Example
%   % Test int.
%   none = []; s = 1; v = 1:10; m = ones(3, 3); type = 'int';
%   fprintf('Test 1: are all elements of the objects %s?  ', type);
%   if istypeof(none, type) && istypeof(s, type) && istypeof(v, type) && istypeof(m, type)
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end
%
%   % Test real.
%   none = []; s = 2.4; v = linspace(1, 10, 17); m = [pi 1; exp(1) 0]; type = 'real';
%   fprintf('Test 2: are all elements of the objects %s?  ', type);
%   if istypeof(none, type) && istypeof(s, type) && istypeof(v, type) && istypeof(m, type)
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end
%
%   % Test Axis.
%   none = []; s = Axis.x; v = Axis.elems; m = [Axis.x Axis.y; Axis.y Axis.z]; type = 'Axis';
%   fprintf('Test 3: are all elements of the objects %s?  ', type);
%   if istypeof(none, type) && istypeof(s, type) && istypeof(v, type) && istypeof(m, type)
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end
%
%   % Test a cell type.
%   none = {}; c = {1, [2 3]; 6, [7 8; 9 10]}; type = 'intcell';
%   fprintf('Test 4: are all objects %s?  ', type);
%   if istypeof(none, type) && istypeof(c, type)
%       fprintf('Yes.\n');
%   else
%       fprintf('No.\n');
%   end

%%% See Also
% <istypesizeof.html istypesizeof>


function truth = istypeof(obj, typename)

chkarg(ischar(typename), '"typename" should be string.');
iscelltype = false;
if ~isempty(strfind(typename, 'cell'))  && isequal(typename(end-3:end), 'cell')  % typename ends with 'cell'
	typename = typename(1:end-4);
	iscelltype = true;
end

if isempty(obj)
	truth = ~xor(iscelltype, iscell(obj));  % true if "typename" is a cell type and obj == {}, or if "typename" is a non-cell type and obj == []
elseif iscelltype
	truth = iscell(obj);
	obj = obj(:);
	for i = 1:length(obj)
		truth = truth && istypeof(obj{i}, typename);
	end
elseif isequal(typename, 'int')
	truth = isnumeric(obj) && isint(obj);  % note isint('abc') == true
elseif isequal(typename, 'real')
	truth = isnumeric(obj) && isreal(obj);  % note isreal('abc') == true
elseif isequal(typename, 'complex')
	truth = isnumeric(obj);
elseif isequal(typename, 'arbitrary')
	truth = true;
else
	truth = isa(obj, typename);
end

function truth = isint(obj)

truth = isreal(obj) && all(mod(obj(:),1)==0);
