function [Table1,h_uifig] = three_fluorophore_FRET_Table1(FP_original_spectra,FP_photophysical_params)

% prepare output table
varTypes = {'string', 'single', 'single', 'single'};
varNames = {'FRET_comb', 'R0', 'r_10', 'r_10_2acc'};
Table1 = table('Size',[5 4],'VariableTypes',varTypes,'VariableNames',varNames);
Table1.FRET_comb = {'mTRQ2-mVEN';'mVEN-mRFP';'mTRQ2-mRFP';...
    'mTRQ2-mVEN-mRFP_(middle)';'mTRQ2-mVEN-mRFP_(random)'};

% photophysical and physicochemical parameters
MAC3 = FP_photophysical_params.epsilon(3);
MAC2 = FP_photophysical_params.epsilon(2);
QY1 = FP_photophysical_params.QY(1);
QY2 = FP_photophysical_params.QY(2);
n = 1.4;

% extract spectra
spec_wave = FP_original_spectra.wavelength;
spec_1_em = FP_original_spectra.mTurquoise2_em;
spec_2_abs = FP_original_spectra.mVenus_abs;
spec_2_em = FP_original_spectra.mVenus_em;
spec_3_abs = FP_original_spectra.mRFP1_abs;

% R0 for 1->2
R0_12 = calc_R0(spec_wave, spec_1_em, spec_2_abs, QY1, MAC2, n, 2/3);

% R0 for 2->3
R0_23 = calc_R0(spec_wave, spec_2_em, spec_3_abs, QY2, MAC3, n, 2/3);

% R0 for 1->3
R0_13 = calc_R0(spec_wave, spec_1_em, spec_3_abs, QY1, MAC3, n, 2/3);

% store in table
Table1.R0(1:3) = [R0_12; R0_23; R0_13];

% distance for each FRET pair corresponding to 10% energy
% transfer efficiency
r_10_percent_12 = R0_12*((1-0.1)/0.1)^(1/6);
r_10_percent_23 = R0_23*((1-0.1)/0.1)^(1/6);
r_10_percent_13 = R0_13*((1-0.1)/0.1)^(1/6);

% store in table
Table1.r_10(1:3) = [r_10_percent_12; r_10_percent_23; r_10_percent_13];

% Second, add the intermediate F2 in a linear arrangement:
% 1 - 2 - 3
% For this very simple estimation, we assume that cross-excitation is
% negligible and 3 can be excited by either 1->3 or 1->2->3. The distance
% between F1 and F3 is increased, while F2 is kept exactly in the middle.

% set the fluorophore distances  
r_13 = 0.1:0.01:20;
r_12 = 0.5 * r_13;
r_23 = 0.5 * r_13;

% pre-allocate
E_middle = zeros(size(r_13));

% loop through all input distances
for ii = 1:numel(r_13)
    E_middle(ii) = calcFRETcascade(r_12(ii), r_13(ii), r_23(ii),...
        R0_12, R0_13, R0_23);
end

% estimate the distance that corresponds to 10% energy transfer efficiency
idx = interp1(E_middle, 1:length(E_middle), 0.1, 'nearest');
r_10_percent_13_middle = r_13(idx);

% store in table
Table1.r_10(4) = r_10_percent_13_middle;

% The intermediate fluorophore will not be equidistant, but most
% likely closer to one of the other fluorophores. In the following, its
% position will be randomized 1000 times for each distance between 1 and 3,
% and then averaged.
E_random = zeros(size(r_13,1), 1000);

for ii = 1:numel(r_13)
    
    for j = 1:1000
        r_12 = 0.1 + (r_13(ii)-0.1).*rand(1);
        r_23 = r_13(ii) - r_12;
        E_random(ii,j) = calcFRETcascade(r_12, r_13(ii), r_23,...
            R0_12, R0_13, R0_23);
    end
    
end

E_random = mean(E_random,2);

% estimate the distance that corresponds to 10% energy transfer efficiency
idx = interp1(E_random, 1:length(E_random), 0.1, 'nearest');
r10_percent_13_random = r_13(idx);

% store in table
Table1.r_10(5) = r10_percent_13_random;

% For membrane-bound fluorophores, this situation is changed. This 2D
% compartment offers less space for proteins, making it more likely that
% one donor is surrounded by multiple acceptors.
% The following calculations and approximations are based on:
% Bunt et al., FRET from single to multiplexed signaling events,
% Biophys. Rev. 9, 2017, 119-129,
% https://doi.org/10.1007/s12551-017-0252-z
%
% The antenna effect: the donor has a higher chance to depopulate the
% excited state via FRET. Although all proteins reside in a two-dimensional
% membrane, we still assume a mostly linear arrangement for the individual
% membrane protein complexes (the calculations would be different, if one
% donor would approach a 2D plane of acceptors with a high surface
% density). Since we only consider that one donor is within range of
% several identical and equidistant acceptors, we can use the approximation
% that R0^6 is effectively multiplied by the number of acceptors, thus
% changing the distance term to (1/n)r^6.
%
% The FRET surplus effect: instead of the average dipole orientation factor
% of 2/3, the donor has a higher probability to find an acceptor that has a
% more favorable relative orientation during each excitation cycle, thus
% raising the average orientation factor kappa_squared. The data provided
% in figure 5 (right panel) of the above publication can be fitted with:
% average_kappa_squared = A1*exp(-N_Acc/t1) + % A2*exp(-na/t2) + y0, with
% A1 = -0.69677, A2 = -0.45809, t1 = 1.00167, t2 = 3.84384, y0 = 1.27663
% and na as the number of acceptors. The function itself is not based
% the underlying physics, it serves merely as a method to interpolate the
% figure data.

% number of acceptors
na = 2;

% R0 for 1->2 (multiple acceptors)
R0_12_multacc = calc_R0_multacc(spec_wave, spec_1_em, spec_2_abs,...
    QY1, MAC2, n, na);

% R0 for 2->3 (multiple acceptors)
R0_23_multacc = calc_R0_multacc(spec_wave, spec_2_em, spec_3_abs,...
    QY2, MAC3, n, na);

% R0 for 1->3 (multiple acceptors)
R0_13_multacc = calc_R0_multacc(spec_wave, spec_1_em, spec_3_abs,...
    QY1, MAC3, n, na);

% use the adjusted distance dependency for multiple acceptors within range
% of one donor in combination with the adjusted R0 (including the adjusted
% kappa_squared values for multiple acceptors)

% r vs. E curves for all fluorophore FRET pairings and their respective
% distances corresponding to 10% energy transfer efficiency
r = 0.1:0.01:20;

E_12_multacc = R0_12_multacc^6./((1/na).*r.^6 + R0_12_multacc^6);
idx = interp1(E_12_multacc,1:length(r),0.1,'nearest');
r_10_percent_12_multacc = r(idx);
E_23_multacc = R0_23_multacc^6./((1/na).*r.^6 + R0_23_multacc^6);
idx = interp1(E_23_multacc,1:length(r),0.1,'nearest');
r_10_percent_23_multacc = r(idx);
E_13_multacc = R0_13_multacc^6./((1/na).*r.^6 + R0_13_multacc^6);
idx = interp1(E_13_multacc,1:length(r),0.1,'nearest');
r_10_percent_13_multacc = r(idx);

% store in table
Table1.r_10_2acc(1:3) =...
    [r_10_percent_12_multacc;r_10_percent_23_multacc;r_10_percent_13_multacc];

% extending these calculations for the FRET cascade
r_13 = r;

% mVEN in the middle
r_12 = 0.5 * r_13;
r_23 = 0.5 * r_13;
E_multacc_middle = zeros(size(r_13));

for ii = 1:numel(r_13)
    E_multacc_middle(ii) = calcFRETcascade_multacc(r_12(ii),r_13(ii),r_23(ii),...
        R0_12_multacc,R0_13_multacc,R0_23_multacc,na);
end

idx = interp1(E_multacc_middle,1:length(r),0.1,'nearest');
r10_percent_multacc_middle = r(idx);

% store in table
Table1.r_10_2acc(4) = r10_percent_multacc_middle;

% mVEN at a random position
E_multacc_random = zeros(size(r_13,1),1000);

for ii = 1:numel(r_13)
    
    for j = 1:1000
        r_12 = 0.1 + (r_13(ii)-0.1).*rand(1);
        r_23 = r_13(ii) - r_12;
        E_multacc_random(ii,j) = calcFRETcascade_multacc(r_12,r_13(ii),r_23,...
            R0_12_multacc,R0_13_multacc,R0_23_multacc,na);
    end
    
end

E_multacc_random = mean(E_multacc_random,2);
idx = interp1(E_multacc_random,1:length(r),0.1,'nearest');
r10_percent_multacc_random = r(idx);

% store in table
Table1.r_10_2acc(5) = r10_percent_multacc_random;

% plot table to figure
close all force;
h_uifig = uifigure('Color','w','Position',[-664 588 401 159]);
uitable(h_uifig,'Data',Table1,...
    'Position',[11 11 381 139],'ColumnWidth',{180,50,73,73});

end