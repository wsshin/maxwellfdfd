classdef PhysC < handle
	% PhysC is the class defining fundamental physical constants.
	
	properties (Constant)
		% Fundamental constants
		c0 = 2.99792458e8;  % in m/s
		mu0 = 4*pi*1e-7;  % in H/m
		eps0 = 1 / (PhysC.c0^2 * PhysC.mu0);  % in F/m
		eta0 = sqrt(PhysC.mu0 /PhysC.eps0);  % in Ohm
		h = 4.135667662e-15;  % in eV·sec
		hbar = PhysC.h / (2*pi);  % in eV·sec
	end
end

