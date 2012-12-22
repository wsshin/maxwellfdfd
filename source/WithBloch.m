classdef WithBloch < handle
	% Abstract class inherited by source classe that can set the Bloch boundary
	% condition.
	
	properties (Abstract, SetAccess = immutable)
		kBloch
	end
end
