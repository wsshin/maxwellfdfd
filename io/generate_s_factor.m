function s_factor_cell = generate_s_factor(omega, grid3d, m, R)
% For example, s_factor_cell{Axis.x, GT.dual} is the s-factor multiplied to Ex
% (or eps_ww).

chkarg(istypesizeof(omega, 'complex'), '"omega" should be complex.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

if nargin < 3  % no m
	m = 4;
end
chkarg(istypesizeof(m, 'real') && m >= 1, '"m" should be real and at least 1.0.');

if nargin < 4  % no R
	lnR = -16;  % R = 1e-7
else
	chkarg(istypesizeof(R, 'real') && R > 0 && R < 1, '"R" should be positive and less than 1.0.');
	lnR = log(R);  % log() is the natural logarithm
end


s_factor_cell = cell(Axis.count, GT.count);
for w = Axis.elems
	Nw = grid3d.N(w);

	lpmls = grid3d.lpml(w,:);  % locations of PML interfaces
	lpml_n = lpmls(Sign.n); lpml_p = lpmls(Sign.p);

	Lpmls = grid3d.Lpml(w,:);  % thicknesses of PML
	Lpml_n = Lpmls(Sign.n); Lpml_p = Lpmls(Sign.p);

	for g = GT.elems
		s_factor = ones(1, Nw);
		l = grid3d.l{w, g};  % length(l) == Nw, rather than Nw+1.
		for i = 1:Nw
			li = l(i);
			if li < lpml_n
				s_factor(i) = calc_s_factor(lpml_n - li, Lpml_n, omega, m, lnR);
% 				if w == Axis.y
% 					s_factor(i) = -s_factor(i);
% 				end
			elseif li > lpml_p
				s_factor(i) = calc_s_factor(li - lpml_p, Lpml_p, omega, m, lnR);
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
