function R0 = calc_R0(spec_wave,spec_em,spec_abs,QY,MAC,n,kappa_squared)
%CALC_R0 calculates the Förster distance
%
%   inputs:
%   |spec_wave|, Nx1 vector containing the wavelengths for the spectra
%   |spec_em|, Nx1 vector containing the intensity data for the donor
%   emission, corresponding to the wavelength scale |spec_wave| and
%   normalized to the range [0 1]
%   |spec_abs|, Nx1 vector containing the intensity data for the acceptor
%   absorption, corresponding to the wavelength scale |spec_wave| and
%   normalized to the range [0 1]
%   |QY|, quantum yield of the donor
%   |MAC|, molar attenuation coefficient of the acceptor (M-1*cm-1)
%   |n|, refractive index
%   |kappa_squared|, dipole-dipole orientation factor
%
%   output:
%   |R0|, Förster distance (nm) for the chosen donor-acceptor pair 
%
%   2020-08-18 Sven zur Oven-Krockhaus

% cumulative normalization factor
F = sum(spec_em);

% spectral overlap integrand
J = (spec_wave(:,1).^4).*MAC.*spec_em.*spec_abs./ F;

% sum up for integral
J = sum(J);

% calculate R0
R0 = 0.02108 * ((kappa_squared * (n^-4) * QY * J) ^ (1/6));

end