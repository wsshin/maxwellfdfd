classdef Domain < Box
	% Domain is a special subclass of Box for specifying a computation domain.
	
	methods
        function this = Domain(bound, dl)
			super_args = {bound, dl};
			this = this@Box(super_args{:});
		end
	end
end

