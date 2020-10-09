function [fitresult, gof] = fit_spec_em(spec_wave, spec_em)
% FIT_SPEC_EM fits a flourescence emission spectrum
%  using a set of four Gauss functions

% fittype and options
ft = fittype('gauss4');
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0 0 0 0 0 0 0];

% get position of maximum wavelength
[~, max_idx] = max(spec_em);
max_wave = spec_wave(max_idx);

% fit start parameters (empirical approximations)
start_waves = max_wave + [0, 11.9, 41.0, 75.6];
start_amps = [0.50, 0.52, 0.29, 0.07];
start_wids = [15.4, 20.7, 31.5, 47.0];
start_vals = [start_amps; start_waves; start_wids];
opts.StartPoint = start_vals(:);

% Fit model to data
[fitresult, gof] = fit(spec_wave, spec_em, ft, opts);

end