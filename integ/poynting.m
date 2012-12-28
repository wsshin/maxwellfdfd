function S_scalar2d = poynting(varargin)

chkarg(nargin == 5 || nargin == 7, 'five or seven arguments are required.')

iarg = 0;
iarg = iarg + 1; polarization = varargin{iarg};
chkarg(istypesizeof(polarization, 'Axis'), '"argument %d should be "polarization" (instance of Axis).', iarg);

if istypesizeof(varargin{2}, 'Scalar3d')
	iarg = iarg + 1; Ep3d = varargin{iarg};
	chkarg(istypesizeof(Ep3d, 'Scalar3d'), '"argument %d should be "Ep3d" (instance of Scalar3d)."', iarg);

	iarg = iarg + 1; Eq3d = varargin{iarg};
	chkarg(istypesizeof(Eq3d, 'Scalar3d'), '"argument %d should be "Eq3d" (instance of Scalar3d)."', iarg);

	iarg = iarg + 1; Hp3d = varargin{iarg};
	chkarg(istypesizeof(Hp3d, 'Scalar3d'), '"argument %d should be "Hp3d" (instance of Scalar3d)."', iarg);

	iarg = iarg + 1; Hq3d = varargin{iarg};
	chkarg(istypesizeof(Hq3d, 'Scalar3d'), '"argument %d should be "Hq3d" (instance of Scalar3d)."', iarg);

	iarg = iarg + 1; normal_axis = varargin{iarg};
	chkarg(istypesizeof(normal_axis, 'Axis'), 'argument %d should be "normal_axis" (instance of Axis).', iarg);

	iarg = iarg + 1; intercept = varargin{iarg};
	chkarg(istypesizeof(intercept, 'real'), 'argument %d should be "intercept" (real).', iarg);
		
	grid3d = Ep3d.grid3d;
	chkarg(isequal(Eq3d.grid3d, grid3d) && isequal(Hp3d.grid3d, grid3d) && isequal(Hq3d.grid3d, grid3d), ... 
		'instances of Scalar3d do not have same grid3d.');
	
	Ep2d = slice_scalar3d(Ep3d, normal_axis, intercept);
	Eq2d = slice_scalar3d(Eq3d, normal_axis, intercept);
	Hp2d = slice_scalar3d(Hp3d, normal_axis, intercept);
	Hq2d = slice_scalar3d(Hq3d, normal_axis, intercept);

	grid2d = Ep2d.grid2d;
else
	iarg = iarg + 1; Ep2d = varargin{iarg};
	chkarg(istypesizeof(Ep2d, 'Scalar2d'), '"argument %d should be "Ep2d" (instance of Scalar2d)."', iarg);

	iarg = iarg + 1; Eq2d = varargin{iarg};
	chkarg(istypesizeof(Eq2d, 'Scalar2d'), '"argument %d should be "Eq2d" (instance of Scalar2d)."', iarg);

	iarg = iarg + 1; Hp2d = varargin{iarg};
	chkarg(istypesizeof(Hp2d, 'Scalar2d'), '"argument %d should be "Hp2d" (instance of Scalar2d)."', iarg);

	iarg = iarg + 1; Hq2d = varargin{iarg};
	chkarg(istypesizeof(Hq2d, 'Scalar2d'), '"argument %d should be "Hq2d" (instance of Scalar2d)."', iarg);
	
	
	grid2d = Ep2d.grid2d;
	intercept = Ep2d.intercept;
	chkarg(Eq2d.intercept==intercept && Hp2d.intercept==intercept && Hq2d.intercept==intercept, ...
		'instances of Scalar2d do not have same intercept.');
	chkarg(isequal(Eq2d.grid2d, grid2d) && isequal(Hp2d.grid2d, grid2d) && isequal(Hq2d.grid2d, grid2d), ... 
		'instances of Scalar2d do not have same grid2d.');
end

array = real(Ep2d.array .* conj(Hq2d.array) - Eq2d.array .* conj(Hp2d.array)) / 2;
osc = Ep2d.osc;
physQ = PhysQ.S;

S_scalar2d = Scalar2d(array, grid2d, osc, physQ, [physQ.symbol, '_', char(polarization)], intercept);
