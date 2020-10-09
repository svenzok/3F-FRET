function output_spec = interpol_spec(input_spec,interpolation_method,steps)
% INTERPOL_SPEC interpolates spectra
%   inputs:
%   |input_spec| as an Nx2 matrix,
%     with columns: wavelength and emission intensity data.
%   |steps| is the desired distance of wavelength units
%   |interpolation_method| can be chosen as 'pchip', 'spline' or 'linear'
%
%   output:
%   |output_spec|, the interpolated spectrum.

% separate x and y data
x = input_spec(:,1);
y = input_spec(:,2);

% determine the x-data range
xmin = ceil(x(1));
xmax = floor(x(end));
xq = xmin:steps:xmax;

% interpolatation according to chosen method
yq = interp1(x,y,xq',interpolation_method);

% normalization to a maximum intensity of 1 (normalization to the interval
% [0 1] is not valid as incomplete spectra may not drop to the baseline)
normy = yq/max(yq);

% output
output_spec = [xq' normy];

end