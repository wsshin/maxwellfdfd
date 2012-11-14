classdef Domain < Box
	% Domain is a special subclass of Box for specifying a computation domain.
	% Unlike other Shape classes, Domain must take "dl_max" as an argument.
	
	methods
        function this = Domain(bound, dl_max)
			super_args = {bound, dl_max};
			this = this@Box(super_args{:});
		end
	end
end

