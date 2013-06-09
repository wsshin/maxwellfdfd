clear all; close all; clear classes; clc;

%% Gold on pp.130-131 in Section 12 of CRC Handbook of Chemistry and Physics, 92nd Ed., W. M. Haynes (Editor), CRC Press
%% Store data in column vectors.
eV = [
0.10
0.20
0.30
0.40
0.50
0.60
0.70
0.80
0.90
1.00
1.20
1.40
1.60
1.80
2.00
2.10
2.20
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
4.10
4.20
4.30
4.40
4.50
4.60
4.70
4.80
4.90
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
27.20
27.60
28.00
28.40
28.80
29.20
29.60
30.00
];

n = [
8.17 
2.13 
0.99 
0.59 
0.39 
0.28 
0.22 
0.18 
0.15 
0.13 
0.10 
0.08 
0.08 
0.09 
0.13 
0.18
0.24 
0.50 
0.82 
1.24 
1.43 
1.46 
1.50 
1.54 
1.54 
1.54 
1.55 
1.56 
1.58 
1.62 
1.64 
1.63 
1.59 
1.55 
1.51 
1.48 
1.45 
1.41 
1.35 
1.30 
1.27 
1.25 
1.23 
1.22 
1.21 
1.21 
1.21 
1.21 
1.22 
1.24 
1.25 
1.27 
1.30 
1.34 
1.36 
1.38 
1.38 
1.35 
1.31 
1.30 
1.30 
1.31 
1.31 
1.30 
1.31 
1.33 
1.36 
1.37 
1.37 
1.36 
1.35 
1.34 
1.34 
1.34 
1.34 
1.35
1.36 
1.38 
1.39 
1.44 
1.45 
1.42 
1.37 
1.33 
1.29 
1.26 
1.24 
1.22 
1.21 
1.20 
1.19 
1.19 
1.19 
1.19 
1.19 
1.20 
1.21 
1.21 
1.18 
1.14 
1.10 
1.05 
1.00 
0.94 
0.89 
0.85 
0.82 
0.80 
0.80 
0.80 
0.80 
0.82 
0.83 
0.84 
0.85 
0.85 
0.86 
0.86 
0.87 
0.88 
0.88 
0.88 
0.87 
0.86
];

k = [82.83 41.73 27.82 20.83 16.61 13.78 11.75 10.21 9.01 8.03 6.54 5.44 4.56 3.82 3.16 2.84 2.54 1.86 1.59 1.54 1.72 1.77 1.79 1.80 1.81 1.80 1.78 1.76 1.73 1.73 1.75 1.79 1.81 1.81 1.79 1.78 1.77 1.76 1.74 1.69 1.64 1.59 1.54 1.49 1.40 1.33 1.27 1.20 1.14 1.09 1.05 1.01 0.97 0.95 0.95 0.96 0.98 0.99 0.96 0.92 0.89 0.88 0.86 0.83 0.81 0.78 0.78 0.79 0.80 0.80 0.80 0.79 0.77 0.76 0.74 0.73 0.72 0.71 0.71 0.73 0.79 0.84 0.86 0.86 0.86 0.84 0.83 0.81 0.79 0.78 0.76 0.75 0.74 0.74 0.73 0.74 0.76 0.80 0.83 0.85 0.87 0.88 0.88 0.86 0.83 0.79 0.75 0.70 0.66 0.62 0.58 0.56 0.54 0.52 0.51 0.50 0.49 0.49 0.48 0.48 0.48 0.48 0.48 0.48];
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


%% Save data.
save(mfilename, 'eV', 'n', 'k');
