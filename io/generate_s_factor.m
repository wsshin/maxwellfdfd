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

	lpml_n = grid3d.lpml(w,Sign.n); lpml_p = grid3d.lpml(w,Sign.p);  % locations of PML interfaces
	Lpml_n = grid3d.Lpml(w,Sign.n); Lpml_p = grid3d.Lpml(w,Sign.p);  % thicknesses of PML
	m_n = m(w,Sign.n); m_p = m(w,Sign.p);
	lnR_n = lnR(w,Sign.n); lnR_p = lnR(w,Sign.p);

	for g = GT.elems
		s_factor = ones(1, Nw);
		l = grid3d.l{w, g};  % length(l) == Nw, rather than Nw+1.
		for i = 1:Nw
			li = l(i);
			if li < lpml_n
				s_factor(i) = calc_s_factor(lpml_n - li, Lpml_n, omega, m_n, lnR_n);
% 				if w == Axis.y
% 					s_factor(i) = -s_factor(i);
% 				end
			elseif li > lpml_p
				s_factor(i) = calc_s_factor(li - lpml_p, Lpml_p, omega, m_p, lnR_p);
% 				if w == Axis.y
% 					s_factor(i) = -s_factor(i);
% 				end
			end
		end
		
		s_factor_cell{w,g} = s_factor;
	end	
end


function s_factor = calc_s_factor(depth, Lpml, omega, m, lnR)
assert(Lpml > 0);

sigma_max = -(m+1) * lnR/(2*Lpml);  % -(m+1) ln(R)/(2 eta Lpml), where eta = 1 in the unit of eta_0
sigma = sigma_max * (depth/Lpml)^m;

kappa_max = 1;
kappa = 1 + (kappa_max-1) * (depth/Lpml)^m;

ma = m;
amax = 0;
a = amax * (1 - depth/Lpml)^ma;

s_factor = kappa + sigma/(a + 1i*omega);  % s = kappa + sigma/(a + i omega eps), where eps = 1 in the unit of eps_0
