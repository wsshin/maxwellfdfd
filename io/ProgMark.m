% ProgMark is a class that marks the progress of operations.
classdef ProgMark < handle
	
	properties (SetAccess = immutable, GetAccess = private)
		ticID_0  % initial tic ID
	end
	
	properties (Access = private)
		t_prev  % elapsed time at previous measurement
	end
	
	methods
		function this = ProgMark()
			this.ticID_0 = tic;
			this.t_prev = 0;
		end
		
		function mark(this, operation_name)
			t_curr = toc(this.ticID_0);  % elapsed time until now
			fprintf('time elapsed: %s sec in total, %s sec for %s\n', ...
				num2str(t_curr), num2str(t_curr - this.t_prev), operation_name);
			this.t_prev = t_curr;
		end
	end
	
end

