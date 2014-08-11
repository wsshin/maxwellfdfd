clear all; close all; clear classes; clc;

%% Tungsten on pp.148 in Section 12 of CRC Handbook of Chemistry and Physics, 92nd Ed., W. M. Haynes (Editor), CRC Press
%% Store data in column vectors.
eV = [
0.10
0.20
0.25
0.30
0.34
0.38
0.42
0.46
0.50
0.54
0.58
0.62
0.66
0.70
0.74
0.78
0.82
0.86
0.90
0.94
0.98
1.00
1.10
1.20
1.30
1.40
1.50
1.60
1.70
1.80
1.90
2.00
2.10
2.20
2.30
2.40
2.50
2.60
2.70
2.80
2.90
3.00
3.10
3.20
3.30
3.40
3.50
3.60
3.70
3.80
3.90
4.00
4.20
4.40
4.60
4.80
5.00
5.20
5.40
5.60
5.80
6.00
6.20
6.40
6.60
6.80
7.00
7.20
7.40
7.60
7.80
8.00
8.20
8.40
8.60
8.80
9.00
9.20
9.40
9.60
9.80
10.00
10.20
10.40
10.60
10.80
11.00
11.20
11.40
11.60
11.80
12.00
12.40
12.80
13.20
13.60
14.00
14.40
14.80
15.20
15.60
16.00
16.40
16.80
17.20
17.60
18.00
18.40
18.80
19.20
19.60
20.00
20.40
20.80
21.20
21.60
22.00
22.40
22.80
23.20
23.60
24.00
24.40
24.80
25.20
25.60
26.00
26.40
26.80
27.00
27.50
28.00
28.50
29.00
29.50
30.00
31.00
32.00
33.00
34.00
35.00
36.00
37.00
38.00
39.00
40.00
];

n = [
14.06 
3.87 
2.56 
1.83 
1.71 
1.86 
1.92 
1.69 
1.40 
1.23 
1.17 
1.28 
1.45 
1.59 
1.83 
2.12 
2.36 
2.92 
3.11 
3.15 
3.15 
3.14 
3.05 
3.00 
3.12 
3.29 
3.48 
3.67 
3.84 
3.82 
3.70 
3.60 
3.54 
3.49 
3.49 
3.45 
3.38 
3.34 
3.31 
3.31 
3.32 
3.35 
3.39 
3.43 
3.45 
3.39 
3.24 
3.13 
3.05 
2.99 
2.96 
2.95 
3.02 
3.13 
3.24 
3.33 
3.40 
3.27
2.92 
2.43 
2.00 
1.70 
1.47 
1.32 
1.21 
1.12 
1.06 
1.01 
0.98 
0.95 
0.93 
0.94 
0.94 
0.96 
0.99 
1.01 
1.01 
1.02 
1.03 
1.05 
1.09 
1.13 
1.19 
1.24 
1.27 
1.29 
1.28 
1.27 
1.25 
1.22 
1.20 
1.16 
1.10 
1.04 
0.98 
0.94 
0.91 
0.90 
0.90 
0.93 
0.97 
0.98 
0.97 
0.94 
0.90 
0.85 
0.80 
0.74 
0.69 
0.64 
0.60 
0.56 
0.54 
0.52 
0.50 
0.50 
0.49 
0.49
0.49 
0.49 
0.48 
0.49 
0.50 
0.51 
0.53 
0.55 
0.57 
0.59 
0.61 
0.62 
0.64 
0.67 
0.69 
0.71 
0.73 
0.75 
0.78 
0.79 
0.82 
0.84 
0.85 
0.85 
0.84 
0.83 
0.81 
0.80
];

k = [54.71 28.30 22.44 18.32 15.71 13.88 12.63 11.59 10.52 9.45 8.44 7.52 6.78 6.13 5.52 5.00 4.61 4.37 4.44 4.43 4.36 4.32 4.04 3.64 3.24 2.96 2.79 2.68 2.79 2.91 2.94 2.89 2.84 2.76 2.72 2.72 2.68 2.62 2.55 2.49 2.45 2.42 2.41 2.45 2.55 2.66 2.70 2.67 2.62 2.56 2.50 2.43 2.33 2.32 2.41 2.57 2.85 3.27 3.58 3.70 3.61 3.42 3.24 3.04 2.87 2.70 2.56 2.43 2.30 2.18 2.06 1.95 1.86 1.76 1.70 1.65 1.60 1.55 1.50 1.44 1.38 1.34 1.33 1.34 1.36 1.39 1.42 1.44 1.46 1.48 1.48 1.48 1.47 1.44 1.40 1.35 1.28 1.23 1.17 1.13 1.12 1.14 1.17 1.19 1.21 1.21 1.20 1.18 1.15 1.11 1.07 1.02 0.97 0.92 0.87 0.82 0.77 0.73 0.69 0.66 0.62 0.57 0.53 0.49 0.46 0.43 0.40 0.38 0.37 0.36 0.34 0.32 0.31 0.30 0.30 0.29 0.29 0.29 0.28 0.29 0.31 0.32 0.33 0.33 0.33 0.33];
k = k.';

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
