%% vis2d
% Visualize a 2D slice of a 3D field solution with objects.

%%% Syntax
%  vis2d(scalar3d, normal_axis, intercept, [opts])
%  vis2d(scalar3d, normal_axis, intercept, obj_array, [opts])

%%% Description
% |vis2d(scalar3d, normal_axis, intercept)| vizualizes a 2D slice of one of the
% x-, y-, z-components of the E- or H-field stored as an instance of
% <Scalar3d.html |Scalar3d|>.
%
% |visall(scalar3d, obj_array)| visualizes the objects in |obj_array| with
% |scalar3d|.  The elements of |obj_array| are instances of <Object.html
% Object>.
%
% The optional argument |opts| is an options structure that controls the
% behavior of visualization.  The fields of the structure are explained below.
% The default values of the fields are shown in brackets |{}|.
%
% * |opts.withgrid|: |true| or |{false}| to show the grid or not.
% * |opts.withinterp|: |{true}| or |false| to interpolate the field or not.
% * |opts.withpml|: |true| or |{false}| to show the PML regions or not.
% * |opts.withabs|: |true| or |{false}| to show the absolute values of the field
% or not
% * |opts.cscale|: positive value multiplied to the color bar range.  Set to
% values less than |1.0| to saturate colors.  The default value is |1.0|.
% * |opts.isopaque|: |true| or |{false}| to make slices opaque or not.
% * |opts.withobj|: |{true}| or |false| to show the objects or not.

%%% Example
%  [E, H, obj_array] = maxwell_run(...);
%  opts.cscale = 5e-3;
%  opts.wihabs = true;
%  visall(E{Axis.y}, obj_array, opts);

function visall(scalar3d, varargin)

narginchk(1, 3)
chkarg(istypesizeof(scalar3d, 'Scalar3d'), '"scalar3d" should be instance of Scalar3d.');

iarg = 2;
ivararg = 1;
obj_array = [];
if iarg <= nargin && ~istypesizeof(varargin{ivararg}, 'struct')
	obj_array = varargin{ivararg};
	chkarg(istypesizeof(obj_array, 'Object', [1 0]), ...
		'argument %d should be "obj_array" (row vector with Object as elements).', iarg);
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
if no_opts || ~isfield(opts, 'isopaque')
	opts.isopaque = false;
end
if no_opts || ~isfield(opts, 'withobj')
	opts.withobj = true;
end

tv = TotalView3d();
tv.scalar3d = scalar3d;
tv.obj_array = obj_array;	
tv.withgrid = opts.withgrid;
tv.withinterp = opts.withinterp;
tv.withpml = opts.withpml;
tv.withabs = opts.withabs;
tv.cscale = opts.cscale;
tv.isopaque = opts.isopaque;
tv.withobj = opts.withobj;

tv.show();
