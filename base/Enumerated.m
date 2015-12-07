%% Enumerated 
% Abstract superclass for enumerated quantities. 

%%% Description
% |Enumerated| is the superclass of classes representing enumerated quantities.
% Each subclass of |Enumerated| allows only instances predefined by
% |enumeration|.  For more details about |enumeration|, see
% <matlab:web(['jar:file://',fullfile(matlabroot,'help','techdoc','help.jar!','matlab_oop','br4imrp.html')],'-helpbrowser')
% Enumerations>.
%
% |Enumerated| supports conversion of its instance into a string by |char()|,
% and into an integer by |int()|.
%
% Each concrete subclass of |Enumerated| allows comparison between its
% instances, i.e., supports relational operators <, <=, ==, >=, and >.  Also,
% the result consistent with this comparison is returned by |sort()| on a vector
% of these instances.

%%% Concrete Subclasses
% * <Axis.html |Axis|>: Catesion axes
% * <BC.html |BC|>: boundary conditions
% * <Dir.html |Dir|>: horizontal and vertical directions
% * <GT.html |GT|>: grid kinds
% * <PhysQ.html |PhysQ|>: physical quantities
% * <PML.html |PML|>: PML kinds
% * <Sign.html |Sign|>: negative and positive signs

%%% Methods (static)
% Below, |Enumerated| must be replaced by a concrete subclass.
%
% * |Enumerated.elems()|: returns a row vector of all instances
% * |Enumerated.count()|: returns the number of all instances

%%% Methods (public)
% Below, |enumerated| is an instance of a concrete subclass of |Enumerated|.
%
% * |char(enumerated)|: name of |enumerated|
% * |int(enumerated)|: integer representation of |enumerated|

%%% Example
% Below, the methods of |Enumerated| are demonstrated using a concrete subclass
% <Axis.html |Axis|>.
%
%   % Test conversion to integers and strings.
%   fprintf('Total # of Axis objects: %d\n', Axis.count);
%   for w = Axis.elems
%       fprintf('Name of Axis object %d: %s\n', int(w), char(w));
%   end
%   
%   % Test sort().
%   unsorted = [Axis.z, Axis.x, Axis.y];
%   fprintf('\nBefore sort: ');
%   for w = unsorted
%       fprintf('%s  ', char(w));
%   end
%
%   fprintf('\nAfter sort: ');
%   sorted = sort(unsorted);
%   for w = sorted
%       fprintf('%s  ', char(w));
%   end
%   fprintf('\n');

classdef Enumerated < handle
	properties (SetAccess = immutable, GetAccess = private)
		name
	end
	
	methods (Abstract, Static)
		elems = elems(ind)
		count = count()
	end

%   (Deleted because it uses Axis rather than Enumerated, and also because
%   loading does not generate warnings without loadobj(); only saveobj() seems
%   necessary.)
% 	methods(Static)
% 		function obj = loadobj(S)
% 			% If saveobj() and loadobj() are not implemented, a warning is
% 			% issued when loading a saved Enumerated object.
% 			obj = Axis(char(S));
% 		end
% 	end		
	
	methods
		function this = Enumerated(name)
			chkarg(ischar(name), '"name" should be string.');
			this.name = name;
		end
		
		function S = saveobj(this)
			% If saveobj() and loadobj() are not implemented, a warning is
			% issued when loading a saved Enumerated object.
			S.name = char(this);
		end
		
		function name = char(this)
			% This function returned name = this.name originally, but it is
			% modified to handle arguments given as arrays.
% 			name = this.name;
			name = cell2mat(arrayfun(@(x) x.name, this, 'UniformOutput', false));
		end
		
		function n = int(this)
			% This function returned n = find(this.elems==this) originally, but
			% it is modified to handle arguments given as arrays.
% 			[~, n] = ismember(this, this.elems, 'legacy');
% 			[~, n] = ismember(this, this.elems);
			n = arrayfun(@(x) find(this(1).elems==x), this);
		end
		
		function ind = subsindex(this)
			ind = int(this) - 1;
		end
		
		function c = plus(this, another)
				c = int(this) + double(another);
		end
				
		function [sorted, ind] = sort(this, varargin)
			n = length(this);
			nums = NaN(1, n);
			for i = 1:n
				nums(i) = int(this(i));
			end
			[~, ind] = sort(nums, varargin{:}); 
			sorted = this(ind);
		end

		function truth = lt(this, another)
			cn = class(this);  % class name
			chkarg(istypesizeof(another, cn), 'cannot compare %s with %s.', cn, class(another));
			truth = int(this) < int(another);
		end
		
		function truth = le(this, another)
			cn = class(this);  % class name
			chkarg(istypesizeof(another, cn), 'cannot compare %s with %s.', cn, class(another));
			truth = int(this) <= int(another);
		end
		
		function truth = gt(this, another)
			cn = class(this);  % class name
			chkarg(istypesizeof(another, cn), 'cannot compare %s with %s.', cn, class(another));
			truth = int(this) > int(another);
		end
		
		function truth = ge(this, another)
			cn = class(this);  % class name
			chkarg(istypesizeof(another, cn), 'cannot compare %s with %s.', cn, class(another));
			truth = int(this) >= int(another);
		end
	end
end
