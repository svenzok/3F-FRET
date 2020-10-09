function [fitresult, gof] = fit_spec_em_logn(spec_wave, spec_em)
% FIT_SPEC_EM_LOGN fits a flourescence emission spectrum
%  using a set of four lognormal functions

% fittype and options
logn_fun_str = [...
    'A1/(w1*x)*exp(-(log(x/xc1))^2/(2*w1^2))+',...
    'A2/(w2*x)*exp(-(log(x/xc2))^2/(2*w2^2))+',...
    'A3/(w3*x)*exp(-(log(x/xc3))^2/(2*w3^2))+',...
    'A4/(w4*x)*exp(-(log(x/xc4))^2/(2*w4^2))'];
    
ft = fittype(logn_fun_str,'independent','x','dependent','y');
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';


% get position of maximum wavelength
[~, max_idx] = max(spec_em);
max_wave = spec_wave(max_idx);

% start parameters (empirical approximations)
start_A = [5.7, 7.0, 6.5, 4.2];
start_w = [0.020, 0.025, 0.038, 0.077];
start_xc = max_wave + [-6.4, 5.8, 35.0, 60.9];
opts.StartPoint = [start_A, start_w, start_xc];
opts.Lower = [1 0.5 2 0.1 0.01 0.01 0.02 0.02 start_xc-40];
opts.Upper = [15 15 15 15 0.05 0.05 0.06 0.10 start_xc+40];

% Fit model to data
[fitresult, gof] = fit(spec_wave, spec_em, ft, opts);

end

