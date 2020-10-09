function output_int = align_spec(input_spec,spec_wave)
% ALIGN_SPEC aligns spectra to a certain wavelength scale
%   inputs:
%   |input_spec| as an Nx2 matrix,
%     with columns: wavelength and emission intensity data.
%   |spec_wave| is the desired wavelength scale
%   
%   |input_spec| and |spec_wave| must have the same wavelength intervals
%   and |input_wave| must fit into the range given by |input_spec|
%   
%   outputs:
%   |output_int|, the aligned emission intensity as an Nx1 vector, with
%     the intensity data shifted to the correct position in relation to
%     the new wavelength scale and missing data replaced with NaN.

% separate x and y data
x = input_spec(:,1);
y = input_spec(:,2);

% pre-allocate the new y data
output_int = NaN(numel(spec_wave),1);

% get range for insertion
[~,ia,ib] = intersect(spec_wave,x);

% insert
output_int(ia) = y(ib);

end