clear all; close all; clear classes; clc;

%% Aluminum on pp.126-127 in Section 12 of CRC Handbook of Chemistry and Physics, 92nd Ed., W. M. Haynes (Editor), CRC Press
%% Store data in column vectors.
eV = [
0.040
0.050
0.060
0.070
0.080
0.090
0.100
0.125
0.150
0.175
0.200
0.250
0.300
0.350
0.400
0.500
0.600
0.700
0.800
0.900
1.000
1.100
1.200
1.300
1.400
1.500
1.600
1.700
1.800
1.900
2.000
2.200
2.400
2.600
2.800
3.000
3.200
3.400
3.600
3.800
4.000
4.200
4.400
4.600
4.800
5.000
6.000
6.500
7.000
7.500
8.000
8.500
9.000
9.500
10.000
10.500
11.000
11.500
12.000
12.500
13.000
13.500
14.000
14.200
14.400
14.600
14.800
15.000
15.200
15.400
15.600
15.800
16.000
16.200
16.400
16.750
17.000
17.250
17.500
17.750
18.000
18.500
19.000
19.500
20.000
20.500
21.000
21.500
22.000
22.500
23.000
23.500
24.000
24.500
25.000
25.500
26.000
27.000
28.000
29.000
30.000
35.000
40.000
45.000
50.000
55.000
60.000
65.000
70.000
72.500
75.000
77.500
80.000
85.000
90.000
95.000
100.000
110.000
120.000
130.000
140.000
150.000
160.000
170.000
180.000
190.000
200.000
220.000
240.000
260.000
280.000
300.000
];

n = [
98.595
74.997
62.852
53.790
45.784
39.651
34.464
24.965
18.572
14.274
11.733
8.586
6.759
5.438
4.454
3.072
2.273
1.770
1.444
1.264
1.212
1.201
1.260
1.468
2.237
2.745
2.625
2.143
1.741
1.488
1.304
1.018
0.826
0.695
0.598
0.523
0.460
0.407
0.363
0.326
0.294
0.267
0.244
0.223
0.205
0.190
0.130
0.110
0.095
0.082
0.072
0.063
0.056
0.049
0.044
0.040
0.036
0.033
0.033
0.034
0.038
0.041
0.048
0.053
0.058
0.067
0.086
0.125
0.178
0.234
0.280
0.318
0.351
0.380
0.407
0.448
0.474
0.498
0.520
0.540
0.558
0.591
0.620
0.646
0.668
0.689
0.707
0.724
0.739
0.753
0.766
0.778
0.789
0.799
0.809
0.817
0.826
0.840 
0.854 
0.865 
0.876 
0.915 
0.940 
0.957 
0.969 
0.979 
0.987 
0.995 
1.006 
1.025 
1.011 
1.008 
1.007 
1.007 
1.005 
0.999 
0.991 
0.994 
0.991 
0.987 
0.989 
0.990 
0.989 
0.989 
0.990 
0.990 
0.991 
0.992 
0.993 
0.993 
0.994 
0.995
];

k = [
203.701 
172.199 
150.799 
135.500 
123.734 
114.102 
105.600 
89.250 
76.960 
66.930 
59.370 
48.235 
40.960 
35.599 
31.485 
25.581 
21.403 
18.328 
15.955 
14.021 
12.464 
11.181 
10.010 
8.949 
8.212 
8.309 
8.597 
8.573 
8.205 
7.821 
7.479
6.846 
6.283 
5.800 
5.385 
5.024 
4.708 
4.426 
4.174 
3.946 
3.740 
3.552 
3.380 
3.222 
3.076 
2.942 
2.391 
2.173 
1.983 
1.814 
1.663 
1.527 
1.402 
1.286 
1.178 
1.076 
0.979 
0.883 
0.791 
0.700 
0.609 
0.517 
0.417 
0.373
0.327 
0.273 
0.211 
0.153 
0.108 
0.184 
0.073 
0.065 
0.060 
0.055 
0.050 
0.045 
0.042 
0.040 
0.038 
0.036 
0.035 
0.032 
0.030 
0.028 
0.027 
0.025 
0.024 
0.023 
0.022 
0.021 
0.021 
0.020 
0.019 
0.018 
0.018 
0.017 
0.016
0.015
0.014
0.014
0.013
0.010
0.008
0.007
0.006
0.005
0.004
0.004
0.004
0.004
0.024
0.025
0.024
0.028
0.031
0.036
0.030
0.025
0.024
0.021
0.016
0.015
0.014
0.011
0.010
0.009
0.007
0.006
0.005
0.004
0.003
0.002
];

%% Convert the photon energies (in eV) to the wavelengths (in nm).
wvlen = PhysC.h * PhysC.c0 * 1e9 ./ eV;

% %% Reverse the data order, and make them column vectors.
% eV = eV(end:-1:1).';
% n = n(end:-1:1).';
% k = k(end:-1:1).';
% wvlen = wvlen(end:-1:1).';

%% Calculate the permittivity from n and k following the exp(+i w t) time dependence.
eps = (n - 1i*k).^2;

%% Plot real(eps) and imag(eps).  Compare with Fig.3 on p.4375 of Johnson.

nk_wvlen = 1;
eps_eV = 2;
eps_wvlen = 3;
eps_omega = 4;
q_omega = 5;
plotstyle = eps_eV;
figure;
switch plotstyle
    case nk_wvlen  % plot n and k
        loglog(wvlen, n, 'o-', wvlen, k, 'o-')
%         plot(wvlen, n, wvlen, k)
        legend('n', 'k', 'Location', 'SouthEast');
        xlabel 'wavelength (nm)'
%         axis([1e2 1e4 1e-2 1e2])
    case eps_eV  % plot real(eps) and -imag(eps)
        loglog(eV, abs(real(eps)), 'o-', eV, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'Photon Energy (eV)'
%         axis([0.5 6.5 -7 7]);
    case eps_wvlen
        plot(wvlen, real(eps), 'o-', wvlen, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'wavelength (nm)'
%         axis([1e2 1e4 1e-2 1e2])
    case eps_omega
        plot(2*pi./wvlen, real(eps), 'o-', 2*pi./wvlen, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel '\omega (c/nm)'
    case q_omega
        omega = 2*pi./wvlen;
        deps1 = real(eps(2:end)) - real(eps(1:end-1));
        domega = omega(2:end) - omega(1:end-1);
        eps1_inter = (real(eps(2:end)) + real(eps(1:end-1)))/2;
        omega_inter = (omega(2:end) + omega(1:end-1))/2;
        numer = eps1_inter + omega_inter .* (deps1./domega);
        denom = -(imag(eps(2:end)) + imag(eps(1:end-1)));  % extra factor 2
        plot(2*pi./omega_inter, numer./denom, 'o-');
        xlabel 'wavelength (nm)'
        ylabel 'electric Q';
end
title(mfilename)

%% Save data.
save(mfilename, 'eV', 'n', 'k');
