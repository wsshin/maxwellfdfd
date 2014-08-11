clear all; close all; clear classes; clc;

%% Silver on p. 147 in Section 12 of CRC Handbook of Chemistry and Physics, 92nd Ed., W. M. Haynes (Editor), CRC Press
%% Store data in column vectors.
eV = [
0.10
0.15
0.20
0.25
0.30
0.35
0.40
0.45
0.50
0.60
0.70
0.80
0.90
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
3.85
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
20.60
21.20
21.60
22.00
22.40
22.80
23.20
23.60
24.00
24.5
25.0
25.5
26.0
26.5
27.0
27.5
28.0
28.5
29.0
30.0
];

n = [
5.03 
3.00 
2.12 
2.05 
6.39 
2.74
2.49 
3.35 
4.43 
4.71 
4.38 
4.04 
3.80 
3.62 
3.47 
3.35 
3.28 
3.17 
2.98 
2.74 
2.54 
2.36 
2.22 
2.11 
2.01 
1.92 
1.86 
1.81 
1.78 
1.75 
1.71 
1.68 
1.63 
1.59 
1.55 
1.50 
1.44 
1.37 
1.30 
1.24 
1.17 
1.11 
1.08 
1.06 
1.04 
1.05 
1.13 
1.17 
1.21 
1.24 
1.27 
1.17 
1.24 
1.21 
1.15 
1.11 
1.08 
1.04 
1.05 
1.06 
1.07 
1.11 
1.09 
1.11 
1.10 
1.10
1.08 
1.04 
1.02 
1.00 
0.97 
0.95 
0.94 
0.91 
0.89 
0.86 
0.85 
0.81 
0.80 
0.79 
0.81 
0.81 
0.79 
0.78 
0.77 
0.76 
0.76 
0.76 
0.77 
0.77 
0.79 
0.79 
0.79 
0.83 
0.84 
0.87 
0.90 
0.93 
0.94 
0.94 
0.95 
0.96 
0.97 
0.98 
0.98 
1.00 
0.99 
0.99 
0.98 
0.98 
0.97 
0.96 
0.95 
0.92 
0.91 
0.91 
0.89 
0.89 
0.88 
0.86 
0.85 
0.84 
0.82 
0.83 
0.84
];

k = [23.38 15.72 11.34 8.10 9.94 6.21 4.68 3.25 3.22 3.77 3.89 3.82 3.65 3.52 3.40 3.30 3.25 3.28 3.32 3.30 3.23 3.11 2.99 2.88 2.77 2.67 2.56 2.47 2.39 2.34 2.29 2.25 2.21 2.17 2.15 2.12 2.09 2.06 2.01 1.96 1.90 1.83 1.78 1.73 1.62 1.45 1.33 1.29 1.23 1.21 1.20 1.16 1.21 1.22 1.21 1.18 1.14 1.06 1.02 0.97 0.95 0.94 0.92 0.93 0.94 0.95 0.95 0.96 0.95 0.94 0.93 0.91 0.90 0.88 0.88 0.85 0.83 0.79 0.76 0.72 0.69 0.69 0.68 0.67 0.65 0.55 0.52 0.48 0.45 0.42 0.38 0.36 0.32 0.31 0.28 0.27 0.25 0.25 0.24 0.23 0.24 0.25 0.25 0.27 0.27 0.29 0.31 0.31 0.32 0.33 0.33 0.34 0.35 0.35 0.34 0.33 0.33 0.33 0.32 0.31 0.30 0.29 0.26 0.25 0.22];
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
