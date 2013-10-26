clear all; close all; clear classes; clc;

%% Palik's GaAs, p.434
% The data here from Ref. [10] measured only k values, so nothing to take from
% this page.

%% Convert the wavelengths to the photon energies.
eVum = PhysC.h * PhysC.c0 * 1e6;  % (photon energy in eV) = (eVum) / (wavelength in microns);


%% Palik's GaAs, p.435 (using Refs. [9,19])
eV1 = [23:-1:12, 10.8, 10:-1:7].';
n1 = [1.037 1.043 1.058 1.025 0.981 0.936 0.889 0.850 0.836 0.840 0.864 0.895 0.923 0.913 0.901 0.899 1.063].';
k1 = [0.228 0.243 0.245 0.212 0.234 0.267 0.324 0.411 0.503 0.602 0.709 0.791 0.881 0.974 1.136 1.435 1.838].';


%% Palik's GaAs, p.436 (continuing to use Refs. [9,19], and then Ref. [8] from 6 eV)
eV2 = [
6.6
6.3
6.2
6.0
5.94
5.90
5.84
5.80
5.74
5.70
5.64
5.60
5.54
5.50
5.44
5.40
5.34
5.30
5.24
5.20
5.14
5.10
5.04
5.00
4.94
4.90
4.84
4.80
4.74
4.70
4.64
4.60
4.54
4.50
4.44
4.40
4.34
4.30
4.24
4.20
4.14
4.10
4.04
4.00
3.94
3.90
];

n2 = [1.247 1.441 1.424 1.264 1.287 1.288 1.304 1.311 1.321 1.325 1.339 1.349 1.367 1.383 1.410 1.430 1.468 1.499 1.552 1.599 1.699 1.802 2.044 2.273 2.654 2.890 3.187 3.342 3.511 3.598 3.708 3.769 3.850 3.913 4.004 4.015 3.981 3.939 3.864 3.810 3.736 3.692 3.634 3.601 3.559 3.538].';
k2 = [2.047 1.988 1.976 2.472 2.523 2.556 2.592 2.625 2.673 2.710 2.772 2.815 2.885 2.936 3.019 3.079 3.181 3.255 3.384 3.484 3.660 3.795 3.992 4.084 4.106 4.047 3.898 3.770 3.574 3.452 3.280 3.169 3.018 2.919 2.715 2.563 2.368 2.260 2.132 2.069 2.001 1.969 1.935 1.920 1.907 1.904].';


%% Palik's GaAs, p.437 (continuing to use Ref. [8])
eV3 = [
3.84
3.80
3.74
3.70
3.64
3.60
3.54
3.50
3.48
3.46
3.44
3.42
3.40
3.38
3.36
3.34
3.32
3.30
3.28
3.26
3.24
3.22
3.20
3.18
3.16
3.14
3.12
3.10
3.08
3.06
3.04
3.02
3.00
2.98
2.96
2.94
2.92
2.90
2.88
2.86
2.84
2.82
2.80
2.78
2.76
2.74
];

n3 = [
3.512
3.501 
3.488 
3.485 
3.489 
3.495 
3.513 
3.531 
3.541 
3.553 
3.566 
3.580 
3.596 
3.614 
3.635 
3.657 
3.681 
3.709 
3.740 
3.776 
3.818 
3.871 
3.938 
4.023 
4.126 
4.229 
4.313 
4.373 
4.413 
4.439 
4.462 
4.483 
4.509 
4.550 
4.626 
4.755 
4.917 
5.052 
5.107 
5.102 
5.065 
5.015 
4.959 
4.902 
4.845 
4.793
];

k3 = [1.905 1.909 1.920 1.931 1.950 1.965 1.992 2.013 2.024 2.036 2.049 2.062 2.076 2.091 2.107 2.123 2.142 2.162 2.183 2.207 2.232 2.260 2.288 2.307 2.304 2.270 2.212 2.146 2.082 2.029 1.988 1.961 1.948 1.952 1.967 1.960 1.885 1.721 1.529 1.353 1.206 1.088 0.991 0.912 0.846 0.789].';


%% Palik's GaAs, p.438 (continuing to use Ref. [8])
eV4 = [
2.72
2.70
2.68
2.66
2.64
2.62
2.60
2.58
2.56
2.54
2.52
2.50
2.48
2.46
2.44
2.42
2.40
2.38
2.36
2.34
2.32
2.30
2.28
2.26
2.24
2.22
2.20
2.18
2.16
2.14
2.12
2.10
2.08
2.06
2.04
2.02
2.00
1.98
1.96
1.94
1.92
1.90
1.88
1.86
1.84
1.82
1.80
];

n4 = [
4.741 
4.694 
4.649 
4.605 
4.567 
4.525 
4.492 
4.456 
4.423 
4.392 
4.362 
4.333 
4.305 
4.279 
4.254 
4.229 
4.205 
4.183 
4.162 
4.141 
4.120 
4.100 
4.082 
4.063 
4.045 
4.029 
4.013 
3.998 
3.983 
3.968 
3.954 
3.940 
3.927 
3.914 
3.902 
3.890 
3.878 
3.867 
3.856 
3.846 
3.836 
3.826 
3.817 
3.809 
3.799 
3.792 
3.785
];

k4 = [0.739 0.696 0.659 0.626 0.595 0.569 0.539 0.517 0.497 0.476 0.458 0.441 0.426 0.411 0.398 0.385 0.371 0.359 0.347 0.337 0.327 0.320 0.308 0.301 0.294 0.285 0.276 0.266 0.257 0.251 0.245 0.240 0.232 0.228 0.223 0.213 0.211 0.203 0.196 0.187 0.183 0.179 0.173 0.173 0.168 0.158 0.151].';


%% Palik's GaAs, p.439 (continuing to use Ref. [8], and then Ref. [5] from 1.40 eV)
% Note that Ref. [5] did not measure k values becaues they are too small.
eV51 = [
1.78
1.76
1.74
1.72
1.70
1.68
1.66
1.64
1.62
1.60
1.58
1.56
1.54
1.52
1.50
];

n51 = [3.779 3.772 3.762 3.752 3.742 3.734 3.725 3.716 3.707 3.700 3.693 3.685 3.679 3.672 3.666].';
k51 = [0.152 0.134 0.127 0.118 0.112 0.105 0.101 0.097 0.093 0.091 0.089 0.087 0.085 0.083 0.080].';

eV52 = [	
1.35
1.30
1.25
1.20
1.15
1.10
1.05
1.00
0.95
0.90
0.85
0.80
0.75
0.70
0.65
];

n52 = [3.5690 3.5388 3.5138 3.4920 3.4724 3.4546 3.4383 3.4232 3.4094 3.3965 3.3847 3.3737 3.3636 3.3543 3.3457].';
k52 = zeros(size(n52));

eV5 = [eV51; eV52];
n5 = [n51; n52];
k5 = [k51; k52];


%% Palik's GaAs, p.440 (continuing to use Ref. [5], and adding k values from Ref. [4] from 0.1 eV)
% Note that Ref. [5] did not measure k values becaues they are too small.
% Adding the k values from Ref. [4] to n values from Ref. [5] may not be
% consistent, but there is no other way around because no references measured
% the n and k values at the same time.  For silicon, I ignored the k values in
% such a case because k values are too small, but silicon does not have
% interband transitions for those frequency range.  For GaAs, on the other hand,
% we cannot ignore the k values, because there is an interband transition
% between 1e1 and 1e2 eV as shown in Fig. 2 on p.433.
eV61 = [
0.60
0.55
0.50
0.45
0.40
0.35
0.30
0.29
0.28
0.27
0.26
0.25
0.24
0.23
0.22
0.21
0.20
0.19
0.18
0.17
0.16
0.15
0.14
0.13
0.12
0.11
];

n61 = [3.3378 3.3306 3.3240 3.3180 3.3125 3.3075 3.3027 3.3017 3.3008 3.2998 3.2988 3.2978 3.2968 3.2954 3.2946 3.2934 3.2921 3.2907 3.2891 3.2874 3.2854 3.2831 3.2803 3.2769 3.2727 3.2671].';
k61 = zeros(size(n61));

eV62 = [0.10 0.09 0.08].';
n62 = [3.2597 3.2493 3.2336].';
k62 = [4.93e-6 1.64e-5 2.83e-5].';

eV6 = [eV61; eV62];
n6 = [n61; n62];
k6 = [k61; k62];


%% Palik's GaAs, p.441, before Refs. [2,3] have both the n and k values
% For n, I continue to use Ref. [5], and then Refs. [2,3] from 0.0645 eV.  For
% k, I continue to use Ref. [4].
eV7 = [0.070 0.065 0.0645 0.0595 0.057 0.054 0.052 0.0495 0.047].';
n7 = [3.2081 3.1886 3.205 3.176 3.157 3.126 3.100 3.058 2.997].';
k7 = 1e-3 * [0.232 8.27 6.27 3.18 4.50 4.56 3.43 2.07 2.33].';


%% Palik's GaAs, pp.441-442 after Refs. [2,3] have both the n and k values
% From this poist, it seems like cm^-1 is the unit Refs. [2,3] used.
invcm81 = [360 350 340 330 320 315].';
n81 = [2.913 2.851 2.770 2.659 2.495 2.380].';
k81 = 1e-2 * [0.650 0.840 1.13 1.59 2.43 3.14].';

invcm82 = [
310 
305 
300 
299 
298 
297 
296 
295 
294 
293 
292 
291 
290 
289 
288 
287 
286 
285 
284 
283 
282 
281 
280 
279 
278 
277 
276 
275 
274 
273 
272 
271 
270 
269.5 
269 
268.5 
268 
267.5 
267 
266.5 
266 
265 
264 
263 
262 
261 
260
];

n82 = [2.229 2.020 1.707 1.623 1.529 1.422 1.298 1.151 0.975 0.761 0.536 0.391 0.323 0.291 0.275 0.267 0.266 0.270 0.278 0.290 0.307 0.329 0.358 0.395 0.443 0.506 0.592 0.713 0.890 1.17 1.64 2.56 4.66 6.63 9.30 11.6 12.4 12.0 11.3 10.6 9.90 8.87 8.12 7.55 7.11 6.76 6.47].';
k82 = [
0.0421 
0.0602 
0.0959 
0.107
0.122
0.141
0.166
0.201
0.257
0.357
0.552
0.826
1.09
1.34
1.57
1.79
2.01
2.23
2.46
2.69
2.94
3.20
3.48
3.79
4.14
4.52
4.97
5.49
6.12
6.91
7.94
9.32
11.0 
11.6 
11.3
9.37 
6.74 
4.66 
3.30 
2.43 
1.87 
1.20 
0.844 
0.629 
0.489 
0.393 
0.323
];

invcm8 = [invcm81; invcm82];
eV8 = eVum ./ ((1./invcm8) .* 1e4);
n8 = [n81; n82];
k8 = [k81; k82];


%% Palik's GaAs, p.443
% I continue to use Refs. [2,3].  After 60 microns, Refs. [2,3] start to use
% microns rather than cm^-1 as units.  After 100 microns, Refs. [2,3] did not
% measure k values, but according to Fig. 2 on p.433 the wavelengths are longer
% than the wavelength for interband transition.  Therefore, I think I can assume
% that GaAs is lossless for these wavelengths.
invcm91 = [
255 
250 
240 
230 
220 
210 
200
];

um91 = (1./invcm91) .* 1e4;
n91 = [5.57 5.08 4.57 4.30 4.13 4.02 3.93].';
k91 = 1e-2 * [15.3 9.02 4.26 2.49 1.63 1.15 0.849].';

um92 = [60 70 80 90 100 150 200 300 500 1000].';
n92 = [3.77 3.71 3.681 3.662 3.650 3.623 3.615 3.611 3.607 3.606].';
k92 = 1e-3 * [3.89 3.79 1.84 0 0 0 0 0 0 0].';  % 3.79e-3 for 70 microns is from Ref. [1]

um9 = [um91; um92];
eV9 = eVum ./ um9;
n9 = [n91; n92];
k9 = [k91; k92];


%% Aggregate data.
eV = [eV1; eV2; eV3; eV4; eV5; eV6; eV7; eV8; eV9];
n = [n1; n2; n3; n4; n5; n6; n7; n8; n9];
k = [k1; k2; k3; k4; k5; k6; k7; k8; k9];


%% Reverse the data order.
eV = eV(end:-1:1);
n = n(end:-1:1);
k = k(end:-1:1);
wvlen = eVum ./ eV;

%% Calculate the permittivity from n and k following the exp(+i w t) time dependence.
eps = (n - 1i*k).^2;

%% Plot n and k.  Compare with Fig.10 on p.753 of Palik.
nk_wvlen = 1;
eps_eV = 2;
eps_wvlen = 3;
abseps_eV = 4;
plotstyle = nk_wvlen;
switch plotstyle
    case nk_wvlen  % plot n and k
        loglog(wvlen, n, '-', wvlen, k, '-')
        %plot(wvlen, n, wvlen, k)
        legend('n', 'k', 'Location', 'SouthEast');
        xlabel 'wavelength (nm)'
        %axis([1e2 1e4 1e-2 1e2])
    case eps_eV  % plot real(eps) and -imag(eps)
        plot(eV, real(eps), 'o-', eV, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'Photon Energy (eV)'
        axis([10 30 -2 1]);
    case eps_wvlen
        plot(wvlen, real(eps), 'o-', wvlen, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'wavelength (nm)'
        %axis([0 1e2 -1e1 1e1])
    case abseps_eV  % plot real(eps) and -imag(eps)
        loglog(eV, abs(real(eps)), 'o-', eV, -imag(eps), 'o-')
        legend('abs(\epsilon_1)', '\epsilon_2', 'Location', 'SouthEast');
		line([10 10], [2e-6 0.5e2], 'color', 'r', 'linestyle', '--');
		line([20 20], [2e-6 0.5e2], 'color', 'r', 'linestyle', '--');
		line([30 30], [2e-6 0.5e2], 'color', 'r', 'linestyle', '--');
        xlabel 'Photon Energy (eV)'
end


%% Save data.
% save(mfilename, 'eV', 'n', 'k');
