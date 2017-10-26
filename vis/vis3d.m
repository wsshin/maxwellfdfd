%% vis3d
% Visualize 2D slices of a field solution with objects and sources in a 3D
% domain.

%%% Syntax
%  vis3d(scalar3d, x_intercepts, y_intercepts, z_intercepts, [opts])
%  vis3d(scalar3d, x_intercepts, y_intercepts, z_intercepts, obj_array, [opts])
%  vis3d(scalar3d, x_intercepts, y_intercepts, z_intercepts, obj_array, src_array, [opts])
%  vis3d(scalar3d, [opts])
%  vis3d(scalar3d, obj_array, [opts])
%  vis3d(scalar3d, obj_array, src_array, [opts])

%%% Description
% |vis3d(scalar3d, x_intercepts, y_intercepts, z_intercepts)| vizualizes
% x-, y-, z-normal slices of a 3D array of data stored in an instance of
% <Scalar3d.html |Scalar3d|>. |x_intercepts|, |y_intercepts|,
% |z_intercepts| are arrays of real numbers indicating the x-, y-,
% z-intercepts of the slices.  Any of the three arrays of intercepts can be
% empty arrays (|[]|) to show no slices normal to the corresponding
% directions. Any of the three arrays of intercepts can be |NaN| to show
% slices normal to the corresponding directions at default locations.
%
% |vis3d(scalar3d)| is equivalent to |vis3d(scalar3d, NaN, NaN, NaN)|, and
% it vizualizes x-, y-, z-normal slices at default locations.  The default
% locations are automatically set at the most prominent locations.
%
% |vis3d(..., obj_array)| visualizes the objects in |obj_array| with the
% slices.  The elements of |obj_array| are instances of <EMObject.html
% |EMObject|>.
%
% |vis3d(..., obj_array, src_array)| visualizes the objects and sources in
% |obj_array| and |src_array| with the slices.  The elements of |src_array|
% are instances of <Source.html |Source|>.
%
% The optional argument |opts| is an options structure that controls the
% behavior of visualization.  The fields of the structure are explained below.
% The default values of the fields are shown in brackets |{}|.
%
% * |opts.withgrid|: |true| or |{false}| to show the grid or not.
% * |opts.withinterp|: |{true}| or |false| to interpolate the field or not.
% * |opts.withpml|: |true| or |{false}| to show the PML regions or not.
% * |opts.withabs|: |true| or |{false}| to show the absolute values of the field
% or not.
% * |opts.cscale|: positive value that controls the color bar range.  Set to
% values less than |1.0| to saturate colors.  The default value is |1.0|.
% * |opts.cmax|: another positive value that controls the color bar range; the
% positive maximum of the color bar range is |cscale * cmax|.  If it is set to
% |NaN|, which is the default, then |cscale * (maximum amplitude of the plotted
% data)| is used for the maximum of the color bar range.
% * |opts.withcolorbar|: |{true}| or |false| to show the colorbar or not.
% * |opts.withobjsrc|: |true| or |false| to show the objects and sources or not.
% The default values is |true| when the number of objects is small (<= 20), but
% |false| when the number of objects is large (> 20).  To show only sources
% without objects, provide an empty |obj_array| in |vis3d()|.
% * |opts.phase|: additional phase angle |phi|.  The field is visualized with an
% additional factor |exp(i*phi)| multiplied.  Not useful if |opts.withabs =
% true|.

%%% Example
%   [E, H, obj_array, src_array] = maxwell_run({ARGUMENTS});
%   opts.cscale = 5e-3;
%   opts.wihabs = true;
%   vis3d(E{Axis.y}, 0, [], 0, obj_array, src_array, opts);

function vis3d(varargin)

iarg = 0;

% Set up an instance of Scalar3d.
iarg = iarg + 1; 
arg = varargin{iarg};
chkarg(iarg <= nargin && (istypesizeof(arg, 'Scalar3d')), ...
		'argument %d should be instance of Scalar3d.', iarg);
scalar3d = arg;

intercept = {NaN, NaN, NaN};
iarg = iarg + 1;
arg = varargin{iarg};
if iarg <= nargin && istypesizeof(arg, 'real', [1 0])
	chkarg(iarg+2 <= nargin, 'x-, y-, z-intercepts are needed.');
	for w = Axis.elems
		chkarg(istypesizeof(arg, 'real', [1 0]), '"argument %d should be "intercept" (real array).', iarg);
		intercept{w} = arg;
		iarg = iarg + 1;
		arg = varargin{iarg};
	end
end

obj_array = [];
if iarg <= nargin && ~istypesizeof(arg, 'struct')
	obj_array = varargin{iarg};
	chkarg(istypesizeof(obj_array, 'EMObject', [1 0]), ...
		'argument %d should be "obj_array" (row vector with EMObject as elements).', iarg);
	iarg = iarg + 1;
	arg = varargin{iarg};
end

src_array = [];
if iarg <= nargin && ~istypesizeof(arg, 'struct')
	src_array = varargin{iarg};
	chkarg(istypesizeof(src_array, 'Source', [1 0]), ...
		'argument %d should be "src_array" (row vector with Source as elements).', iarg);
	iarg = iarg + 1;
	arg = varargin{iarg};
end

no_opts = true;
if iarg <= nargin
	opts = arg;
	chkarg(istypesizeof(opts, 'struct'), 'argument %d should be "opts" (struct).', iarg);
	no_opts = false;
	iarg = iarg + 1;
end

chkarg(iarg > nargin, 'Too many arguments.');


if no_opts || ~isfield(opts, 'withgrid')
	opts.withgrid = false;
end
if no_opts || ~isfield(opts, 'withinterp')
	opts.withinterp = true;
end
if no_opts || ~isfield(opts, 'withpml')
	opts.withpml = false;
end
if no_opts || ~isfield(opts, 'withabs')
	opts.withabs = false;
end
if no_opts || ~isfield(opts, 'cscale')
	opts.cscale = 1.0;
end
if no_opts || ~isfield(opts, 'cmax')
	opts.cmax = NaN;
end
if no_opts || ~isfield(opts, 'withcolorbar')
	opts.withcolorbar = true;
end
if no_opts || ~isfield(opts, 'withobjsrc')
	if length(obj_array) <= 20
		opts.withobjsrc = true;
	else
		opts.withobjsrc = false;
		warning('Maxwell:vis', ['for fast plotting, object and sources are not drawn by default if they are more than 20.\n', ...
			'Suggestion: set "opts.withobjsrc = true" to draw objects.']);
	end
end
if no_opts || ~isfield(opts, 'phase')
	opts.phase = 0;
end


p = Painter3d();
p.scalar3d = scalar3d;
p.obj_array = obj_array;	
p.src_array = src_array;
p.withgrid = opts.withgrid;
p.withinterp = opts.withinterp;
p.withpml = opts.withpml;
p.withabs = opts.withabs;
p.cscale = opts.cscale;
p.cmax = opts.cmax;
p.withcolorbar = opts.withcolorbar;
p.phase_angle = opts.phase;

p.init_display();
for w = Axis.elems
	if isnan(intercept{w})
		p.draw_slice(gca, w);
	else
		for i = intercept{w}
			p.draw_slice(gca, w, i);
		end
	end
end

if opts.withobjsrc
	p.draw_objsrc();
end
