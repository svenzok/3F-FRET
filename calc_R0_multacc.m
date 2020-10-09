function R0_multacc = calc_R0_multacc(spec_wave,spec_em,spec_abs,QY,MAC,n,na)
%CALC_R0_MULTACC calculates R0 in case of multi-acceptor conditions
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
%   |na|, number of identical, equidistant acceptors within the Förster
%   distance of the donor
%
%   output:
%   |R0_multacc|, Förster distance (nm) for the chosen
%   donor-acceptor pair and number of acceptors
%
% The approximation for kappa_squared for multi-acceptor conditions is
% based on:
% Bunt et al., FRET from single to multiplexed signaling events,
% Biophys. Rev. 9, 2017, 119-129,
% https://doi.org/10.1007/s12551-017-0252-z
%
% The FRET surplus effect: instead of the average dipole orientation factor
% of 2/3, the donor has a higher probability to find an acceptor that has a
% more favorable relative orientation during each excitation cycle, thus
% raising the average orientation factor kappa_squared. The data provided
% in in figure 5 (right panel) of the above mentioned publication can be
% fitted with the function: average_kappa_squared = A1*exp(-N_Acc/t1) +
% A2*exp(-na/t2) + y0, with
% A1 = -0.69677, A2 = -0.45809, t1 = 1.00167, t2 = 3.84384, y0 = 1.27663
% and |na| as the number of acceptors. The function itself is not based
% the underlying physics, it serves merely as a method to interpolate the
% figure data. For this approximation, viable acceptors have to be inside
% the Förster radius of the donor.
%
% 2020-08-18 Sven zur Oven-Krockhaus

% cumulative normalization factor
F = sum(spec_em);

% spectral overlap integrand
J = (spec_wave(:,1).^4).*MAC.*spec_em.*spec_abs./ F;

% sum up for integral
J = sum(J);

% approximation for plot in Fig.5 (right panel) in Bunt et al. 2017
A1 = -0.69677; A2 = -0.45809; t1 = 1.00167; t2 = 3.84384; y0 = 1.27663;
average_kappa = A1*exp(-na/t1) + A2*exp(-na/t2) + y0;

% calculate R0 (multiple acceptors within the donor's R0)
R0_multacc = 0.02108*((average_kappa*(n^-4)*QY*J)^(1/6));

end