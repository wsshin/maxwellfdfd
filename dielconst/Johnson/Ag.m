clear all; close all; clear classes; clc;

%% Ag in Table I on p.4374 of Johnson and Christy, 1972 Phys. Rev. B
eV = [0.64 0.77 0.89 1.02 1.14 1.26 1.39 1.51 1.64 1.76 1.88 2.01 2.13 2.26 2.38 2.50 2.63 2.75 2.88 3.00 3.12 3.25 3.37 3.50 3.62 3.74 3.87 3.99 4.12 4.24 4.36 4.49 4.61 4.74 4.86 4.98 5.11 5.23 5.36 5.48 5.60 5.73 5.85 5.98 6.10 6.22 6.35 6.47 6.60];
n = [0.24 0.15 0.13 0.09 0.04 0.04 0.04 0.04 0.03 0.04 0.05 0.06 0.05 0.06 0.05 0.05 0.05 0.04 0.04 0.05 0.05 0.05 0.07 0.10 0.14 0.17 0.81 1.13 1.34 1.39 1.41 1.41 1.38 1.35 1.33 1.31 1.30 1.28 1.28 1.26 1.25 1.22 1.20 1.18 1.15 1.14 1.12 1.10 1.07];
k = [14.08 11.85 10.10 8.828 7.795 6.992 6.312 5.727 5.242 4.838 4.483 4.152 3.858 3.586 3.324 3.093 2.869 2.657 2.462 2.275 2.070 1.864 1.657 1.419 1.142 0.829 0.392 0.616 0.964 1.161 1.264 1.331 1.372 1.387 1.393 1.389 1.378 1.367 1.357 1.344 1.342 1.336 1.325 1.312 1.296 1.277 1.255 1.232 1.212];

%% Convert the photon energies to the wavelengths.
wvlen = PhysC.h * PhysC.c0 * 1e9 ./ eV;

%% Make the row vectors into column vectors to be used in interp1q().
eV = eV.';
n = n.';
k = k.';
wvlen = wvlen.';

%% Calculate the permittivity from n and k following the exp(+i w t) time dependence.
eps = (n - 1i*k).^2;

%% Plot real(eps) and imag(eps).  Compare with Fig.3 on p.4375 of Johnson.

nk_wvlen = 1;
eps_eV = 2;
eps_wvlen = 3;
eps_omega = 4;
q_omega = 5;
plotstyle = eps_eV;
switch plotstyle
    case nk_wvlen  % plot n and k
        loglog(wvlen, n, 'o-', wvlen, k, 'o-')
        %plot(wvlen, n, wvlen, k)
        legend('n', 'k', 'Location', 'SouthEast');
        xlabel 'wavelength (nm)'
        %axis([1e2 1e4 1e-2 1e2])
    case eps_eV  % plot real(eps) and -imag(eps)
        plot(eV, real(eps), 'o-', eV, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'Photon Energy (eV)'
%         axis([0.5 6.5 -7 7]);
    case eps_wvlen
        plot(wvlen, real(eps), 'o-', wvlen, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'wavelength (nm)'
        %axis([1e2 1e4 1e-2 1e2])
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
