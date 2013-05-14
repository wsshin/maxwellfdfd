clear all; close all; clear classes; clc;

%% Ag in Table I on p.4374 of Johnson and Christy, 1972 Phys. Rev. B
eV = [0.64 0.77 0.89 1.02 1.14 1.26 1.39 1.51 1.64 1.76 1.88 2.01 2.13 2.26 2.38 2.50 2.63 2.75 2.88 3.00 3.12 3.25 3.37 3.50 3.62 3.74 3.87 3.99 4.12 4.24 4.36 4.49 4.61 4.74 4.86 4.98 5.11 5.23 5.36 5.48 5.60 5.73 5.85 5.98 6.10 6.22 6.35 6.47 6.60];
n = [0.92 0.56 0.43 0.35 0.27 0.22 0.17 0.16 0.14 0.13 0.14 0.21 0.29 0.43 0.62 1.04 1.31 1.38 1.45 1.46 1.47 1.46 1.48 1.50 1.48 1.48 1.54 1.53 1.53 1.49 1.47 1.43 1.38 1.35 1.33 1.33 1.32 1.32 1.30 1.31 1.30 1.30 1.30 1.30 1.33 1.33 1.34 1.32 1.28];
k = [13.78 11.21 9.519 8.145 7.150 6.350 5.663 5.083 4.542 4.103 3.697 3.272 2.863 2.455 2.081 1.833 1.849 1.914 1.948 1.958 1.952 1.933 1.895 1.866 1.871 1.883 1.898 1.893 1.889 1.878 1.869 1.847 1.803 1.749 1.688 1.631 1.577 1.536 1.497 1.460 1.427 1.387 1.350 1.304 1.277 1.251 1.226 1.203 1.188];

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
        axis([0.5 6.5 -7 7]);
    case eps_wvlen
        plot(wvlen, real(eps), 'o-', wvlen, -imag(eps), 'o-')
        legend('\epsilon_1', '\epsilon_2', 'Location', 'SouthEast');
        xlabel 'wavelength (nm)'
        %axis([1e2 1e4 1e-2 1e2])
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
