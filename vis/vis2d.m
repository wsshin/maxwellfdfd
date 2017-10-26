%% vis2d
% Visualize a 2D slice of a 3D field solution with objects and sources.

%%% Syntax
%  vis2d(scalar3d, normal_axis, intercept, [opts])
%  vis2d(scalar3d, normal_axis, intercept, obj_array, [opts])
%  vis2d(scalar3d, normal_axis, intercept, obj_array, src_array, [opts])
%  vis2d(scalar2d, [opts])
%  vis2d(scalar2d, obj_array, [opts])
%  vis2d(scalar2d, obj_array, src_array, [opts])

%%% Description
% |vis2d(scalar3d, normal_axis, intercept)| vizualizes a 2D slice of one of the
% x-, y-, z-components of the E- or H-field stored as an instance of
% <Scalar3d.html |Scalar3d|>.  |normal_axis| is an instance of <Axis.html Axis>
% the represents the axis normal to the slice, and |intercept| is the intercept
% of the slice on the normal axis.
%
% |vis2d(scalar2d)| vizualizes a 2D slice stored as an instance of
% <Scalar2d.html |Scalar2d|>.
%
% |vis2d(..., obj_array)| visualizes the objects in |obj_array| with the slice
% of a field.  The elements of |obj_array| are instances of <EMObject.html
% |EMObject|>.
%
% |vis2d(..., obj_array, src_array)| visualizes the objects and sources in
% |obj_array| and |src_array| with the slice of a field.  The elements of
% |src_array| are instances of <Source.html |Source|>.
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
% without objects, provide an empty |obj_array| in |vis2d()|.
% * |opts.phase|: additional phase angle |phi|.  The field is visualized with an
% additional factor |exp(i*phi)| multiplied.  Not useful if |opts.withabs =
% true|.

%%% Example
%   [E, H, obj_array, src_array] = maxwell_run({ARGUMENTS});
%   opts.cscale = 5e-3;
%   opts.wihabs = true;
%   vis2d(E{Axis.y}, Axis.z, 0, obj_array, src_array, opts);

function vis2d(varargin)

iarg = 0;

% Set up an instance of Scalar2d.
iarg = iarg + 1; 
chkarg(iarg <= nargin && (istypesizeof(varargin{iarg}, 'Scalar2d') || istypesizeof(varargin{iarg}, 'Scalar3d')), ...
		'argument %d should be instance of either Scalar2d or Scalar3d.', iarg);
arg = varargin{iarg};

if istypesizeof(arg, 'Scalar2d')
	scalar2d = arg;
else
	scalar3d = arg;
	iarg = iarg + 1;
	chkarg(iarg <= nargin && istypesizeof(varargin{iarg}, 'Axis'), '"argument %d should be "normal_axis" (instance of Axis).', iarg);
	normal_axis = varargin{iarg};

	iarg = iarg + 1; 
	chkarg(iarg <= nargin && istypesizeof(varargin{iarg}, 'real'), '"argument %d should be "intercept" (real).', iarg);
	intercept = varargin{iarg};

	scalar2d = slice_scalar3d(scalar3d, normal_axis, intercept);
end

iarg = iarg + 1;
obj_array = [];
if iarg <= nargin && ~istypesizeof(varargin{iarg}, 'struct')
	obj_array = varargin{iarg};
	chkarg(istypesizeof(obj_array, 'EMObject', [1 0]), ...
		'argument %d should be "obj_array" (row vector with EMObject as elements).', iarg);
	iarg = iarg + 1;
end

src_array = [];
if iarg <= nargin && ~istypesizeof(varargin{iarg}, 'struct')
	src_array = varargin{iarg};
	chkarg(istypesizeof(src_array, 'Source', [1 0]), ...
		'argument %d should be "src_array" (row vector with Source as elements).', iarg);
	iarg = iarg + 1;
end

no_opts = true;
if iarg <= nargin
	opts = varargin{iarg};
	chkarg(istypesizeof(opts, 'struct'), 'argument %d should be "opts" (struct).', iarg);
	no_opts = false;
	iarg = iarg + 1;
end

chkarg(iarg > nargin, 'Too many arguments.');


if no_opts || ~isfield(opts, 'isswapped')
	opts.isswapped = false;
end
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
if no_opts || ~isfield(opts, 'linewidth')
	opts.linewidth = 1.0;
end


p = Painter2d();
p.scalar2d = scalar2d;
p.obj_array = obj_array;	
p.src_array = src_array;
p.isswapped = opts.isswapped;
p.withgrid = opts.withgrid;
p.withinterp = opts.withinterp;
p.withpml = opts.withpml;
p.withabs = opts.withabs;
p.cscale = opts.cscale;
p.cmax = opts.cmax;
p.withcolorbar = opts.withcolorbar;
p.phase_angle = opts.phase;
p.linewidth = opts.linewidth;

p.init_display();
p.draw_slice();
if opts.withobjsrc
	p.draw_objsrc();
end
