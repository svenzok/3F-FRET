function [ ratios ] = unmix_tripleFRET(spec_input,spec_comp1,spec_comp2,spec_comp3)
% UNMIX_TRIPLEFRET estimates the fluorophore ratios in a measured spectrum
%  by simple linear unmixing
%
%  inputs:
%  |spec_input| as the spectral data for linear unmixing as an Nx2 matrix,
%    with columns: wavelength and emission intensity data.
%  |spec_comp1|, |spec_comp2| and |spec_comp3| as the emission spectra of
%    the single components, each as an Nx2 matrix, with columns:
%    wavelength and emission intensity data.
%  output:
%  |ratios| as an 1x8 vector, with ratios(peak intensities), ratios(areas),
%  area(rawdata), area(fitspectrum)
%
%  The spectra must be measured in the same wavelength units, but do not
%  need to have the same wavelength scale or intervals as they will be
%  interpolated
%
%  The single component spectra should be measured in the same system as
%  the input data (e.g., by overexpression of each fluorescently tagged
%  protein in separate experiments). Using spectra from literature is, in
%  principle, possible, but might lead to less precise results, as a
%  different physicochemical environment can result in spectral shifts or
%  different spectral shapes.
%
%  All spectra should feature a robust signal-to-noise ratio, as this
%  simple unmixing procedure is not suitable for a high, inhomogeneous
%  background signal. If possible, all spectra should cover the full
%  emission range of the respective fluorophores.

% Some workarounds: For this kind of linear unmixing (solving a system of
% linear equations) no missing values are allowed. However, in many cases
% the wavelength ranges of the components and the measured spectrum are not
% exactly congruent. Since the emission spectra typically settle quickly to
% zero on both sides of the main peak, setting the missing values at the
% far ends to zero is often a valid approximation. We use a slightly
% different approach and try to fit the component spectra with a set of
% lognormal functions. Depending on the quality of the fit, the missing
% values are then filled with the simulated ones or (as a fall-back) with
% zeros. This works reasonably well even for truncated spectra, but, of
% course, aquiring complete spectra should be the first priority.

% put the component spectra in a cell array for easier indexing
spec_comp = cell(3,1);
spec_comp{1} = spec_comp1;
spec_comp{2} = spec_comp2;
spec_comp{3} = spec_comp3;

% set the resolution of the wavelength scale
steps = 0.5;

% the measured spectrum defines the wavelength range
[wave_start,wave_end] = bounds(spec_input(:,1));

% define new wavelength scale
new_wave = (ceil(wave_start):steps:floor(wave_end))';

% interpolate the input spectrum to satisfy the chosen step size
spec_input = interpol_spec(spec_input,'pchip',steps);

max_waves = zeros(1,3);
for ii = 1:3
    % interpolate the component spectra to satisfy the chosen step size
    spec_comp{ii} = interpol_spec(spec_comp{ii},'pchip',steps);
    
    % try to fit them
    % (use fit_spec_em as an alternative)
    [fitresult, gof] =...
        fit_spec_em_logn(spec_comp{ii}(:,1),spec_comp{ii}(:,2));
        
    % Use a large wavelength range
    % here to cover the entire visible emission spectrum.
    allvis_wave = (250:steps:800)';
    
    % put the original component spectrum into this range
    spec_comp_allvis = align_spec(spec_comp{ii}, allvis_wave);
    
    % If the fit is satisfactory, replace the missing values with the
    % approximated ones, otherwise with zeros.  
    
    tf = isnan(spec_comp_allvis);
    if gof.adjrsquare >= 0.999
        approx_spec = feval(fitresult,allvis_wave);
        spec_comp_allvis(tf) = approx_spec(tf);
    else
        spec_comp_allvis(tf) = 0;
    end
    
    % area normalize
    spec_comp_allvis = spec_comp_allvis./...
        trapz(allvis_wave,spec_comp_allvis);
        
    % shift to the correct position in the wavelength scale of the measured
    % spectrum and store in the |spec_comp| cell array
    spec_comp{ii} = align_spec([allvis_wave,spec_comp_allvis], new_wave);
    
    % get the emission maximum positions
    [~, max_idx] = max(spec_comp{ii});
    max_waves(ii) = new_wave(max_idx); 
end

% collect the components in one matrix and also add a background vector (constant value)
comps = [ones(size(new_wave)), cell2mat(spec_comp')];

% align the input spectrum, but replace the NaN values with zeros
input_intensity = align_spec(spec_input, new_wave);
input_intensity(isnan(input_intensity)) = 0;

% create weight vectors (regions with less fluorescence, i.e., with a worse
% signal-to-noise ratio, are considered less in the fit).
% To reduce outliers, the input spectrum will be smoothed significantly for
% this operation.
weight = mat2gray(sgolayfilt(input_intensity,7,151));

% unmixing (the weights are multiplied with each component)
warning('off','MATLAB:rankDeficientMatrix');
c = lsqnonneg(weight .* comps, weight .* input_intensity);

% to refine the results, the fit spectrum is defined as a new weight
weight = mat2gray(comps * c);
c = lsqnonneg(weight .* comps, weight .* input_intensity);

%%%%
% insert iteration loop for batch fitting
%%%%

% component ratios and intensities
% -----------------------------------------------------------

    % via the maximum intensities
    
    % only look in a +/-5 wavelength region around the expected maxima
    % positions of the fluorophores
    wave_range = repmat((-5:steps:5)',1,3) + max_waves;
    [~,idx] = ismember(wave_range,new_wave);
    Imax = max(input_intensity(idx));
    
    % calculate the ratios
    ratio_int = Imax/sum(Imax);
    
    % overall intensity of the measured spectrum
    int_total_1 = trapz(spec_input(:,1),spec_input(:,2));

    % via unmixing (preferred method)
    ratio_comp = c(2:4)/sum(c(2:4));

    % overall intensity of the fit spectrum
    int_total_2 = trapz(new_wave,comps*c);


% results vector
% ------------------------------------------------------------
    ratios = [ratio_int, ratio_comp', int_total_1, int_total_2];

% if input was only one measured spectrum, then show plot
if size(ratios,1) == 1
    close all;
    h_fig = figure('Color','w','Position',[-1537, 420, 1217, 464]);
    h_ax = axes(h_fig);
    hold on;
    for a = 1:3
        area(h_ax,new_wave,comps(:,a+1)*c(a+1));
        plot(h_ax,new_wave,comps(:,a+1)*c(a+1));
    end
    plot(h_ax,spec_input(:,1),spec_input(:,2));
    plot(h_ax,new_wave,comps*c);
    
    % format the figure
    area_h = flipud(findobj(h_ax,'Type','area'));
    line_h = flipud(findobj(h_ax,'Type','line'));
    
    EdgeColor_list = num2cell([0 1 1; 1 1 0; 1 0 0], 2);
    FaceColor_list = EdgeColor_list;
        
    set(area_h,{'FaceColor'},FaceColor_list);
    set(line_h(1:3),{'Color'},EdgeColor_list);
    set(area_h,{'EdgeColor','FaceAlpha'},{'none',0.2});
    set(line_h(4), {'LineStyle','Marker','MarkerSize','MarkerEdgeColor'},...
    {'none','o',10,[0.5 0.5 0.5]});
    marker_idx = get_markerindices(spec_input(:,1),spec_input(:,2),...
        80,h_fig.Position(3)/h_fig.Position(4));
    line_h(4).MarkerIndices = marker_idx;
    line_h(5).LineWidth = 2;
    line_h(5).Color = [0 0 0];
    
h_ax.TickDir = 'out';
h_ax.Box = 'off';
h_ax.XLabel.String = 'wavelength [nm]';
h_ax.YLabel.String = 'emission intensity [a.u.]';
h_ax.FontSize = 12;
h_ax.XLabel.FontSize = 15;
h_ax.YLabel.FontSize = 15;
h_ax.TickLength = [0.005 0.04];
h_ax.XLim = spec_input([1 end],1);
h_ax.YLim = [0 1.1].*max(comps*c);
    
end
    
warning('on','MATLAB:rankDeficientMatrix');

end
