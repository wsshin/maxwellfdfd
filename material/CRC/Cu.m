clear all; close all; clear classes; clc;

%% Copper on p.129 in Section 12 of CRC Handbook of Chemistry and Physics, 92nd Ed., W. M. Haynes (Editor), CRC Press
%% Store data in column vectors.
eV = [
0.10
0.50
1.00
1.50
1.70
1.75
1.80
1.85
1.90
2.00
2.10
2.20
2.30
2.40
2.60
2.80
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
6.50
7.00
7.50
8.00
8.50
9.00
9.50
10.00
11.00
12.00
13.00
14.00
14.50
15.00
15.50
16.00
17.00
18.00
19.00
20.00
21.00
22.00
23.00
24.00
25.00
26.00
27.00
28.00
29.00
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
41.00
42.00
43.00
44.00
45.00
46.00
47.00
48.00
49.00
50.00
51.00
52.00
53.00
54.00
55.00
56.00
57.00
58.00
59.00
60.00
61.00
62.00
63.00
64.00
65.00
66.00
67.00
68.00
69.00
70.00
75.00
80.00
85.00
90.00
];

n = [
29.69 
1.71 
0.44 
0.26 
0.22 
0.21 
0.21 
0.22 
0.21 
0.27 
0.47 
0.83 
1.04 
1.12 
1.15 
1.17 
1.18 
1.23 
1.27 
1.31 
1.34 
1.34 
1.42 
1.49 
1.52 
1.53 
1.47 
1.38 
1.28
1.18 
1.10 
1.04 
0.96 
0.97 
1.00 
1.03 
1.03 
1.03 
1.03 
1.04 
1.07 
1.09 
1.08 
1.06 
1.03 
1.01 
0.98 
0.95 
0.91 
0.89 
0.88 
0.88 
0.90 
0.92 
0.94 
0.96 
0.96 
0.92 
0.88 
0.86 
0.85 
0.86 
0.88 
0.89 
0.90 
0.91 
0.92 
0.92 
0.92 
0.93 
0.93 
0.93 
0.94 
0.94 
0.94 
0.95 
0.95 
0.95 
0.95 
0.95 
0.95 
0.95 
0.95 
0.95 
0.96 
0.96 
0.96 
0.96 
0.96
0.96 
0.97 
0.97 
0.97 
0.97 
0.96 
0.96 
0.97 
0.97 
0.97 
0.97 
0.97 
0.97 
0.98 
0.98 
0.97 
0.96
];

k = [71.57 17.63 8.48 5.26 4.43 4.25 4.04 3.85 3.67 3.24 2.81 2.60 2.59 2.60 2.50 2.36 2.21 2.07 1.95 1.87 1.81 1.72 1.64 1.64 1.67 1.71 1.78 1.80 1.78 1.74 1.67 1.59 1.37 1.20 1.09 1.03 0.98 0.92 0.87 0.82 0.75 0.73 0.72 0.72 0.72 0.71 0.69 0.67 0.62 0.56 0.51 0.45 0.41 0.38 0.37 0.37 0.40 0.40 0.38 0.35 0.30 0.26 0.24 0.22 0.21 0.20 0.20 0.19 0.19 0.18 0.17 0.17 0.16 0.16 0.15 0.15 0.15 0.15 0.14 0.14 0.14 0.13 0.13 0.13 0.12 0.12 0.12 0.11 0.11 0.11 0.11 0.11 0.11 0.11 0.10 0.10 0.10 0.10 0.09 0.09 0.09 0.09 0.09 0.09 0.09 0.08];
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
