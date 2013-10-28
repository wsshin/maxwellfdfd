classdef TotalView3d < handle
	% TotalView3d visualizes x-, y-, z-cross sections of Scalar3d.
	
	properties (SetAccess = immutable)
		painter3d
		painter2d  % [painter2d for x-normal slice, y-normal slice, z-normal slice]
	end
	
	properties (Access = private)
		intercept  % [x intercept, y intercept, z intercept]
		hs3d  % [hs3d_x, hs3d_y, hs3d_z]: figure handles of slices intercepting x-, y-, z-axes
		ha2d  % [ha2d_x, ha2d_y, ha2d_z]: axes handles for 2D slices
		ha3d  %  axes handle for 3D view
		hedit  % [hedit_x, hedit_y, hedit_z]: handles for text edit UIs
		hslider  % [hslider_x, hslider_y, hslider_z]: handles for slider UIs
	end
	
	properties (Dependent)
		scalar3d
		obj_array
		src_array
		withgrid
		withinterp
		withpml
		phase_angle
		withabs
		cscale
		cmax
		isopaque
	end
	
	properties
		withobjsrc
	end
	
	methods
		function this = TotalView3d(scalar3d)
			this.painter3d = Painter3d();
			this.painter2d = [Painter2d(), Painter2d(), Painter2d()];
			this.intercept = NaN(1, Axis.count);
			this.hs3d = NaN(1, Axis.count);
			this.ha2d = NaN(1, Axis.count);
			this.hedit = NaN(1, Axis.count);
			this.hslider = NaN(1, Axis.count);
			this.withobjsrc = true;
			
			if nargin >= 1  % scalar3d
				this.scalar3d = scalar3d;
				this.show();
			end
		end
		
		function scalar3d = get.scalar3d(this)
			scalar3d = this.painter3d.scalar3d;
		end
		
		function set.scalar3d(this, scalar3d)
			this.painter3d.scalar3d = scalar3d;
			this.intercept = scalar3d.grid3d.center;
			
			for w = Axis.elems
				warning('off', 'Maxwell:interp');
				this.painter2d(w).scalar2d = slice_scalar3d(scalar3d, w, this.intercept(w));
				warning('on', 'Maxwell:interp');
			end
		end
		
		function obj_array = get.obj_array(this)
			obj_array = this.painter3d.obj_array;
		end
		
		function set.obj_array(this, obj_array)
			this.painter3d.obj_array = obj_array;

			for w = Axis.elems
				this.painter2d(w).obj_array = obj_array;
			end
		end
		
		function src_array = get.src_array(this)
			src_array = this.painter3d.src_array;
		end
		
		function set.src_array(this, src_array)
			this.painter3d.src_array = src_array;

			for w = Axis.elems
				this.painter2d(w).src_array = src_array;
			end
		end
		
		function withgrid = get.withgrid(this)
			withgrid = this.painter3d.withgrid;
		end
		
		function set.withgrid(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.painter3d.withgrid = truth;
			for w = Axis.elems
				this.painter2d(w).withgrid = truth;
			end
		end
		
		function withinterp = get.withinterp(this)
			withinterp = this.painter3d.withinterp;
		end
		
		function set.withinterp(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.painter3d.withinterp = truth;
			for w = Axis.elems
				this.painter2d(w).withinterp = truth;
			end
		end
		
		function withpml = get.withpml(this)
			withpml = this.painter3d.withpml;
		end
		
		function set.withpml(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.painter3d.withpml = truth;
			for w = Axis.elems
				this.painter2d(w).withpml = truth;
			end
		end
		
		function phase_angle = get.phase_angle(this)
			phase_angle = this.painter3d.phase_angle;
		end
		
		function set.phase_angle(this, phase_angle)
			chkarg(istypesizeof(phase_angle, 'real'), '"phase_angle" should be real.');
			this.painter3d.phase_angle = phase_angle;
			for w = Axis.elems
				this.painter2d(w).phase_angle = phase_angle;
			end
		end
		
		function withabs = get.withabs(this)
			withabs = this.painter3d.withabs;
		end
		
		function set.withabs(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.painter3d.withabs = truth;
			for w = Axis.elems
				this.painter2d(w).withabs = truth;
			end
		end
		
		function cscale = get.cscale(this)
			cscale = this.painter3d.cscale;
		end
		
		function set.cscale(this, cscale)
			chkarg(istypesizeof(cscale, 'real') && cscale > 0, '"cscale" should be positive.');
			this.painter3d.cscale = cscale;
			for w = Axis.elems
				this.painter2d(w).cscale = cscale;
			end
		end

		function cmax = get.cmax(this)
			cmax = this.painter3d.cmax;
		end
		
		function set.cmax(this, cmax)
			chkarg(istypesizeof(cmax, 'real') && (cmax >= 0 || isnan(cmax)), '"cmax" should be positive or NaN.');
			this.painter3d.cmax = cmax;
			for w = Axis.elems
				this.painter2d(w).cmax = cmax;
			end
		end
		
		function isopaque = get.isopaque(this)
			isopaque = this.painter3d.isopaque;
		end
		
		function set.isopaque(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.painter3d.isopaque = truth;
		end
				
		function set.withobjsrc(this, truth)
			chkarg(istypesizeof(truth, 'logical'), '"truth" should be logical.');
			this.withobjsrc = truth;
		end
		
		function show(this, figure_handle)
			if nargin < 3  % no figure_handle
				figure_handle = gcf;
			end
			
			this.init();
			clf(figure_handle);
			set(figure_handle, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
			bcolor = get(figure_handle, 'Color');  % background color
			
			% Set the figure renderer.  MATLAB uses 'zbuffer' by default, but it
			% uses 'opengl' when a figure has alpha data (for transparency)
			% because 'opengl' is the only renderer that supports alpha.  The
			% problem is that 'zbuffer' produces a better figure than 'opengl'
			% when alpha is not used.  TotalView3d needs alpha only in the 3D
			% view panel (or axes), but unfortunately different panels in a
			% single figure cannot have different renderers.  Therefore users
			% may get confused if their 2D slices look very different depending
			% on "totalview.isopacity" value.  To prevent such confusion, I may
			% set the figure renderer always as 'opengl'.
% 			set(figure_handle, 'Renderer', 'opengl');
			
			p3d = this.painter3d;
			p2d = this.painter2d;
			
			% Figure out the longest axis.
			L = p3d.databound(:,Sign.p) - p3d.databound(:,Sign.n);
			[~, i_l] = max(L);  % index for logest axis
			axis_l = Axis.elems(i_l);  % longest axis
			[axis_m, axis_n] = cycle(axis_l);
			
			% Caculate the sizes of panes.
			wm = 0.07;  % margin in width
			wg = 0.07;  % gap in width
			hm = 0.07;  % margin in height
			hg = 0.10;  % gap in height
			w3 = 0.23;  % width of 3rd pane
			h4 = 0.2;  % height of control pane
						
			wr = 1 - 2*wg - 2*wm - w3;  % rest width
			hr = 1 - hg - 2*hm;  % rest height
			
			% [w1 w2]
			w1 = L(axis_l) / (L(axis_l) + L(axis_n)) * wr;
			w2 = L(axis_n) / (L(axis_l) + L(axis_n)) * wr;
			
			% [h2; h1]
			h1 = L(axis_m) / (L(axis_m) + L(axis_n)) * hr;
			h2 = L(axis_n) / (L(axis_m) + L(axis_n)) * hr;
			
			% For each 2D slice view, assign its horizontal and vertical axes.
			p2d(axis_l).isswapped = true;  % [axis_h axis_v] == [axis_n axis_m]
			p2d(axis_m).isswapped = true;  % [axis_h axis_v] == [axis_l axis_n]
			p2d(axis_n).isswapped = false;  % [axis_h axis_v] == [axis_l axis_m]
			
			this.ha2d(axis_l) = axes('Units', 'normalized', 'Position', [(wm+w1+wg) hm w2 h1]);
			this.ha2d(axis_n) = axes('Units', 'normalized', 'Position', [wm hm w1 h1]);
			this.ha2d(axis_m) = axes('Units', 'normalized', 'Position', [wm (hm+h1+hg) w1 h2]);
			
			this.ha3d = axes('Units', 'normalized', 'Position', [(wm+w1+wg) (hm+h1+hg) w2 h2]);

			% Prepare data to show.
			p3d.prep_data();
			this.cmax = p3d.cmax;
			
			% Draw slices in 2D.
			for w = Axis.elems
				p2d(w).init_display(this.ha2d(w));
				p2d(w).draw_slice(this.ha2d(w));
				if this.withobjsrc
					p2d(w).draw_objsrc(this.ha2d(w));
				end
			end
			
			% Draw slices in 3D.
			p3d.init_display(this.ha3d);
% 			upvec = zeros(1, Axis.count);
% 			upvec(axis_n) = 1;
% 			camup(ha3, upvec);
% 			view(ha3, -38.5, 16)

			% Draw a color bar.
 			hc = colorbar('Position', [(wm+w1+wg+w2+wg) (hm+h1+hg) 0.02 h2]);
			format_colorbar(hc, this.scalar3d);

			if this.withobjsrc
				p3d.draw_objsrc();
			end
			
			for w = Axis.elems
				this.hs3d(w) = p3d.draw_slice(this.ha3d, w, this.intercept(w));
			end
			
			% Add sliders and text edits.
			h5 = h4/3;
			w5 = w3 * 0.6;
			w6 = w3 * 0.05;
			h6 = h5 * 0.5;
			for w = Axis.elems
				s = int(w);  % step
				uicontrol('Style', 'text', ...
					'Units', 'normalized', 'Position', [(1-wm-w3) (hm+h1-s*h5) w6 h5], ...
					'String', upper(char(w)), 'HorizontalAlignment', 'left', 'FontSize', 15, ...
					'BackgroundColor', bcolor);
				this.hedit(w) = uicontrol('Style', 'edit', ...
					'Units', 'normalized', 'Position', [(1-wm-w3+w6) (hm+h1-(s-1)*h5-h6+0.08*h5) (w3-w5-w6)*0.8 h6], ...
					'String', num2str(this.intercept(w)), 'HorizontalAlignment', 'left', 'FontSize', 15, ...
					'BackgroundColor', bcolor, ...
					'Callback', @(src,event) edit_updated(this, w));
				this.hslider(w) = uicontrol('Style', 'slider', ...
					'Units', 'normalized', 'Position', [(1-wm-w5) (hm+h1-s*h5) w5 h5], ...
					'Min', p3d.databound(w,Sign.n), 'Max', p3d.databound(w,Sign.p), 'Value', this.intercept(w), ...
					'BackgroundColor', bcolor, ...
					'Callback', @(src,event) slider_updated(this, w));
			end
			
			% Resize the entire figure window slightly.  Without this command
			% the contents of the figure are shown shifted.  Resizing the window
			% recovers the desired layout of the contents.  Note that
			% refresh(figure_handle) does not give the same effect.
			set(figure_handle, 'Units', 'normalized', 'OuterPosition', [0 0 1-eps 1]);
		end
		
		function edit_updated(this, axis)
			delete(this.hs3d(axis));
			cla(this.ha2d(axis));
			h = this.hedit(axis);
			val = str2double(get(h, 'String'));
			if isnan(val)
				val = this.intercept(axis);
				set(h, 'String', num2str(val));
			else
				vmin = this.painter3d.databound(axis, Sign.n);
				vmax = this.painter3d.databound(axis, Sign.p);
				if val < vmin
					val = vmin;
					set(h, 'String', num2str(val));
				elseif val > vmax
					val = vmax;
					set(h, 'String', num2str(val));
				end
				this.intercept(axis) = val;
			end
			set(this.hslider(axis), 'Value', val);
			this.hs3d(axis) = this.painter3d.draw_slice(this.ha3d, axis, val);
			warning('off', 'Maxwell:interp');
			this.painter2d(axis).scalar2d = slice_scalar3d(this.scalar3d, axis, val);
			warning('on', 'Maxwell:interp');
			this.painter2d(axis).draw_slice(this.ha2d(axis));
			if this.withobjsrc
				this.painter2d(axis).draw_objsrc(this.ha2d(axis));
			end
		end
		
		function slider_updated(this, axis)
			delete(this.hs3d(axis));
			cla(this.ha2d(axis));
			val = get(this.hslider(axis), 'Value');
			this.intercept(axis) = val;
			set(this.hedit(axis), 'String', num2str(val));
			this.hs3d(axis) = this.painter3d.draw_slice(this.ha3d, axis, val);
			warning('off', 'Maxwell:interp');
			this.painter2d(axis).scalar2d = slice_scalar3d(this.scalar3d, axis, val);
			warning('on', 'Maxwell:interp');
			this.painter2d(axis).draw_slice(this.ha2d(axis));
			if this.withobjsrc
				this.painter2d(axis).draw_objsrc(this.ha2d(axis));
			end
		end
	end
	
	methods (Access = private)
		function init(this)
			for w = Axis.elems
				h = this.hs3d(w);
				if ~isempty(h) && ishandle(h)
					delete(h);
				end
				
				h = this.ha2d(w);
				if ~isempty(h) && ishandle(h)
					delete(h);
				end
				
				h = this.hedit(w);
				if ~isempty(h) && ishandle(h)
					delete(h);
				end
				
				h = this.hslider(w);
				if ~isempty(h) && ishandle(h)
					delete(h);
				end
			end
				
			h = this.ha3d;
			if ~isempty(h) && ishandle(h)
				delete(h);
			end
		end
	end
end

