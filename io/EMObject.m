classdef EMObject < handle
	% EMObject is the combination of Shape and Material.
	
	properties (SetAccess = immutable)
		shape
		material
	end
	
	methods
		function this = EMObject(shape, material)
			if nargin ~= 0
				chkarg(istypesizeof(shape, 'Shape', [1 0]), '"shape" should be row vector with Shape as elements.');
				chkarg(istypesizeof(material, 'Material'), ...
					'"material" should be instance of Material.');

				n = length(shape);
				this(n) = EMObject;
				for i = 1:n
					this(i).material = material;
					this(i).shape = shape(i);
				end
			end
		end
		
		% A function assigning eps and mu could have been a member function of
		% EMObject, but there are too many subpixel smoothing schemes for eps and
		% mu.  Suppoting all those schemes would make EMObject very heavyweight. A
		% decision of taking out those schemes from EMObject was made to leave
		% EMObject lightweight.
	end
end
