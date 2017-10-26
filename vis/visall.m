%% visall
% Visualize x-, y-, z-normal slices of a 3D field solution with objects and
% sources.

%%% Syntax
%  visall(scalar3d, [opts])
%  visall(scalar3d, obj_array, [opts])
%  visall(scalar3d, obj_array, src_array, [opts])

%%% Description
% |visall(scalar3d)| vizualizes one of the x-, y-, z-components of the E- or
% H-field stored as an instance of <Scalar3d.html |Scalar3d|>.  The field is
% visualized on x-, y-, z-normal slices, which are shown in 3D as well as in 2D.
%
% |visall(..., obj_array)| visualizes the objects in |obj_array| with the slices
% of a field. The elements of |obj_array| are instances of <EMObject.html
% |EMObject|>.
%
% |vis2d(..., obj_array, src_array)| visualizes the objects and sources in
% |obj_array| and |src_array| with the slices of a field.  The elements of
% |src_array| are instances of <Source.html |Source|>.
%
% The optional argument |opts| is an options structure that controls the
% behavior of visualization.  The allowed fields of |opts| are explained below.
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
% * |opts.isopaque|: |{true}| or |false| to show opaque slices in 3D or not.
% * |opts.withobjsrc|: |true| or |false| to show the objects and sources or not.
% The default values is |true| when the number of objects is small (<= 20), but
% |false| when the number of objects is large (> 20).  To show only sources
% without objects, provide an empty |obj_array| in |visall()|.
% * |opts.phase|: additional phase angle |phi|.  The field is visualized with an
% additional factor |exp(i*phi)| multiplied.  Not useful if |opts.withabs =
% true|.

%%% Note
% When |opts.isopaque = false| is used, the 3D view shows the slices nicely with
% some transparency.  However, the 2D views show worse-looking images because
% the figure renderer that supports transparency interpolates colors differently
% from the normal renderer.

%%% Example
%   [E, H, obj_array, src_array] = maxwell_run({ARGUMENTS});
%   opts.cscale = 1e-1;
%   opts.wihabs = true;
%   visall(E{Axis.x}, obj_array, src_array, opts);


function visall(scalar3d, varargin)

narginchk(1, 4)
chkarg(istypesizeof(scalar3d, 'Scalar3d'), '"scalar3d" should be instance of Scalar3d.');

iarg = 2;
ivararg = 1;
obj_array = [];
if iarg <= nargin && ~istypesizeof(varargin{ivararg}, 'struct')
	obj_array = varargin{ivararg};
	chkarg(istypesizeof(obj_array, 'EMObject', [1 0]), ...
		'argument %d should be "obj_array" (row vector with EMObject as elements).', iarg);
	iarg = iarg + 1;
	ivararg = ivararg + 1;
end

src_array = [];
if iarg <= nargin && ~istypesizeof(varargin{ivararg}, 'struct')
	src_array = varargin{ivararg};
	chkarg(istypesizeof(src_array, 'Source', [1 0]), ...
		'argument %d should be "src_array" (row vector with Source as elements).', iarg);
	iarg = iarg + 1;
	ivararg = ivararg + 1;
end

no_opts = true;
if iarg <= nargin
	opts = varargin{ivararg};
	chkarg(istypesizeof(opts, 'struct'), 'argument %d should be "opts" (struct).', iarg);
	no_opts = false;
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
if no_opts || ~isfield(opts, 'isopaque')
	opts.isopaque = true;
end
if no_opts || ~isfield(opts, 'withobjsrc')
	opts.withobjsrc = true;
end
if no_opts || ~isfield(opts, 'phase')
	opts.phase = 0;
end

tv = TotalView3d();
tv.scalar3d = scalar3d;
tv.obj_array = obj_array;
tv.src_array = src_array;
tv.withgrid = opts.withgrid;
tv.withinterp = opts.withinterp;
tv.withpml = opts.withpml;
tv.withabs = opts.withabs;
tv.cscale = opts.cscale;
tv.cmax = opts.cmax;
tv.isopaque = opts.isopaque;
tv.withobjsrc = opts.withobjsrc;
tv.phase_angle = opts.phase;

tv.show();
