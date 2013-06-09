clear all; close all; clear classes; clc;

%% Platinum on pp.138-139 in Section 12 of CRC Handbook of Chemistry and Physics, 92nd Ed., W. M. Haynes (Editor), CRC Press
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
0.55
0.60
0.65
0.70
0.75
0.80
0.85
0.90
0.95
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
3.20
3.40
3.60
3.80
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
16.50
17.00
17.50
18.00
18.50
19.00
19.50
20.00
20.50
21.00
21.50
22.00
22.50
23.00
23.50
24.00
24.50
25.00
25.50
26.00
26.50
27.00
27.50
28.00
28.50
29.00
29.50
30.00
];

n = [
13.21 
8.18 
5.90 
4.70 
3.92 
3.28 
2.81 
3.03 
3.91 
4.58 
5.13 
5.52 
5.71 
5.57 
5.31 
5.05 
4.77 
4.50 
4.25 
3.86 
3.55 
3.29 
3.10 
2.92 
2.76
2.63 
2.51 
2.38 
2.30 
2.23 
2.17 
2.10 
2.03 
1.96 
1.91 
1.87 
1.83 
1.79 
1.75 
1.68 
1.63 
1.58 
1.53 
1.49 
1.45 
1.43 
1.39 
1.38 
1.36 
1.36 
1.36 
1.36 
1.36 
1.38 
1.39 
1.42 
1.45 
1.48 
1.50 
1.50 
1.49 
1.48 
1.48 
1.47 
1.47 
1.47 
1.47 
1.47 
1.48 
1.49 
1.49 
1.49 
1.48 
1.46 
1.43 
1.40 
1.37 
1.35 
1.33 
1.31 
1.30 
1.29 
1.29 
1.29 
1.29
1.29 
1.31 
1.31 
1.31 
1.30 
1.27 
1.27 
1.25 
1.24 
1.24 
1.25 
1.27 
1.31 
1.30 
1.28 
1.23 
1.18 
1.11 
1.03 
0.94 
0.87 
0.81 
0.77 
0.75 
0.74 
0.73 
0.73 
0.73 
0.74 
0.74 
0.74 
0.74 
0.75 
0.75 
0.75 
0.74 
0.73
];

k = [44.72 31.16 23.95 19.40 16.16 13.66 11.38 9.31 7.71 7.14 6.75 6.66 6.83 7.02 7.04 6.98 6.91 6.77 6.62 6.24 5.92 5.61 5.32 5.07 4.84 4.64 4.43 4.26 4.07 3.92 3.77 3.67 3.54 3.42 3.30 3.20 3.10 3.01 2.92 2.76 2.62 2.48 2.37 2.25 2.14 2.04 1.95 1.85 1.76 1.67 1.61 1.54 1.47 1.40 1.35 1.29 1.26 1.24 1.24 1.25 1.23 1.22 1.20 1.18 1.17 1.15 1.14 1.13 1.12 1.11 1.12 1.13 1.15 1.15 1.16 1.15 1.14 1.12 1.10 1.08 1.06 1.04 1.01 1.00 0.97 0.94 0.93 0.93 0.93 0.93 0.93 0.93 0.92 0.89 0.87 0.86 0.85 0.88 0.94 0.99 1.03 1.06 1.09 1.10 1.08 1.04 0.98 0.92 0.87 0.82 0.77 0.73 0.70 0.67 0.65 0.63 0.62 0.60 0.59 0.58 0.58 0.58];
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
