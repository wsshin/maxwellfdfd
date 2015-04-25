% clear all; close all; clear classes; clc;

%% Palik I's Si3N4, p.774
eV = [24 23 22 21 20 19 18 17 16 15 14 13 12 11 10.5 10 9.5 9 8.5 8 7.75 7.5 7.25 7 6.75 6.5 6.25 6 5.75 5.5 5.25 5 4.75 4.5 4.25 4 3.5 3 2.5 2 1.5 1];

n = [0.655 0.625 0.611 0.617 0.635 0.676 0.735 0.810 0.902 1.001 1.111 1.247 1.417 1.657 1.827 2.000 2.162 2.326 2.492 2.651 2.711 2.753 2.766 2.752 2.724 2.682 2.620 2.541 2.464 2.393 2.331 2.278 2.234 2.198 2.167 2.141 2.099 2.066 2.041 2.022 2.008 1.998];

k = [0.420 0.481 0.560 0.647 0.743 0.841 0.936 1.03 1.11 1.18 1.26 1.35 1.43 1.52 1.53 1.49 1.44 1.32 1.16 0.962 0.866 0.750 0.612 0.493 0.380 0.273 0.174 0.102];
k = [k 1e-2 * [5.7 2.9 1.1]];
k = [k 1e-3 * [4.9 1.2]];
k = [k 1e-4 * [2.3]];
k = [k zeros(1,8)];


%% Convert the wavelengths to the photon energies.
wvlen = PhysC.h * PhysC.c0 * 1e6 ./ eV;

%% Reverse the data order.
eV = eV(end:-1:1);
n = n(end:-1:1);
k = k(end:-1:1);
wvlen = wvlen(end:-1:1);

%% Calculate the permittivity from n and k following the exp(+i w t) time dependence.
eps = (n - 1i*k).^2;

%% Plot n and k.  Compare with Fig. 1 on p.1066 of Palik II.
nk_wvlen = 1;
eps_eV = 2;
eps_wvlen = 3;
abseps_eV = 4;
n_wvlen_fit = 5;
plotstyle = n_wvlen_fit;
switch plotstyle
    case nk_wvlen  % plot n and k
        loglog(wvlen, n, '-o', wvlen, k, '-o')
        %plot(wvlen, n, wvlen, k)
        legend('n', 'k', 'Location', 'SouthEast');
        xlabel 'wavelength (µm)'
        %axis([1e2 1e4 1e-2 1e2])
    case eps_eV  % plot real(eps) and -imag(eps)
        plot(eV, real(eps), 'o-', eV, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'Photon Energy (eV)'
        %axis([0.5 6.5 -7 7]);
    case eps_wvlen
        plot(wvlen, real(eps), 'o-', wvlen, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'wavelength (µm)'
        %axis([0 1e2 -1e1 1e1])
    case abseps_eV  % plot real(eps) and -imag(eps)
        loglog(eV, abs(real(eps)), 'o-', eV, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'Photon Energy (eV)'
    case n_wvlen_fit  % plot n and fitting curve taken from http://refractiveindex.info/?shelf=main&book=Si3N4&page=Philipp
		n_fit = sqrt(1 + (2.8939 * wvlen.^2) ./ (wvlen.^2 - 0.13967^2));
        plot(wvlen, n, '-o', wvlen, n_fit, '-')
        legend('n', 'fit', 'Location', 'SouthEast');
        xlabel 'wavelength (µm)'
end
title(mfilename)

%% Save data.
% save(mfilename, 'eV', 'n', 'k');
