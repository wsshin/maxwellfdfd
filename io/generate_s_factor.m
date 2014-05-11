function s_factor_cell = generate_s_factor(omega, grid3d, m, R)
% For example, s_factor_cell{Axis.x, GT.dual} is the s-factor multiplied to Ex
% (or eps_ww).

chkarg(istypesizeof(omega, 'complex'), '"omega" should be complex.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

if nargin < 3  % no m
	m = 4;
end
chkarg(istypeof(m, 'real') && all(m(:) >= 1), 'element of "m" should be real and at least 1.0.');
chkarg(isexpandable2mat(m, Axis.count, Sign.count), ...
	'"m" should be scalar, length-%d vector, or %d-by-%d matrix.', Axis.count, Axis.count, Sign.count);
m = expand2mat(m, Axis.count, Sign.count);

if nargin < 4  % no R
	R = exp(-16);  % R = exp(-16) ~= 1e-7
end
chkarg(istypeof(R, 'real') && all(R(:) > 0) && all(R(:) <= 1), ...
	'element of "R" should be real and between 0.0 and 1.0.');
chkarg(isexpandable2mat(R, Axis.count, Sign.count), ...
	'"R" should be scalar, length-%d vector, or %d-by-%d matrix.', Axis.count, Axis.count, Sign.count);
R = expand2mat(R, Axis.count, Sign.count);
lnR = log(R);  % log() is the natural logarithm

s_factor_cell = cell(Axis.count, GT.count);
for w = Axis.elems
	Nw = grid3d.N(w);

	lpml = grid3d.lpml(w,:);  % locations of PML interfaces
	Lpml = grid3d.Lpml(w,:);  % thicknesses of PML

	for g = GT.elems
		l = grid3d.l{w, g};  % length(l) == Nw, rather than Nw+1.
		ind_pml = {l < lpml(Sign.n), l > lpml(Sign.p)};  % indices of locatiotns indise PML

		s_factor = ones(1, Nw);
		for s = Sign.elems
			s_factor(ind_pml{s}) = calc_s_factor(omega, abs(lpml(s) - l(ind_pml{s})), Lpml(s), m(w,s), lnR(w,s));
		end
				
		s_factor_cell{w,g} = s_factor;
	end	
end


function s_factor = calc_s_factor(omega, depth, Lpml, m, lnR)
% assert(Lpml > 0);

sigma_max = -(m+1) * lnR/(2*Lpml);  % -(m+1) ln(R)/(2 eta Lpml), where eta = 1 in the unit of eta_0
sigma = sigma_max * (depth/Lpml).^m;

kappa_max = 1;
kappa = 1 + (kappa_max-1) * (depth/Lpml).^m;

ma = m;
amax = 0;
a = amax * (1 - depth/Lpml).^ma;

s_factor = kappa + sigma ./ (a + 1i*omega);  % s = kappa + sigma/(a + i omega eps), where eps = 1 in the unit of eps_0
