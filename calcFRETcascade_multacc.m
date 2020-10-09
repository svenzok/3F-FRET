function E = calcFRETcascade_multacc(r_12,r_13,r_23,R0_12_multacc,R0_13_multacc,R0_23_multacc,na)
%CALCFRETCASCADE_MULTACC calculates cascading FRET efficiencies with
%multiple acceptors per donor
%
%   inputs:
%   |r_12|, distance between F1 and F2 (nm)
%   |r_23|, distance between F2 and F3 (nm)
%   |r_13|, distance between F1 and F3 (nm)
%   |R0_12_multacc|, Förster distance between F1 and F2 (nm) 
%   |R0_23_multacc|, Förster distance between F2 and F3 (nm)
%   |R0_13_multacc|, Förster distance between F1 and F3 (nm)
%   |na|, number of identical, equidistant acceptors within the Förster
%   distance of each donor
%
%   all R0 input values must have been calculated considering the adjusted
%   kappa_squared values for multiple acceptors
%
%   output:
%   |E|, energy transfer efficiency from F1 to F3, by summing up the
%   possible pathways of either 1->3 (direct) or 1->2->3 (via the
%   intermediate F2) in case of multiple acceptors

% The calculations for the energy transfer efficiencies are based on:
% Sun et al., Three-Color Spectral FRET Microscopy Localizes Three
% Interacting Proteins in Living Cells,
% Biophys. J. 2010, 99, 4, 1274-1283
% https://doi.org/10.1016/j.bpj.2010.06.004
%
% Their equations were then further adjusted, using the the approximation
% for multi-acceptor conditions in:
% Bunt et al., FRET from single to multiplexed signaling events,
% Biophys. Rev. 9, 2017, 119-129,
% https://doi.org/10.1007/s12551-017-0252-z

f = 1/na;
E_12_multacc = R0_12_multacc.^6*f.*r_13.^6/(R0_12_multacc.^6*f.*r_13^6+R0_13_multacc.^6.*f.*r_12^6+f.*r_12.^6.*f.*r_13^6);
E_13_multacc = R0_13_multacc.^6*f.*r_12.^6/(R0_12_multacc.^6*f.*r_13^6+R0_13_multacc.^6.*f.*r_12^6+f.*r_12.^6.*f.*r_13^6);
E_23_multacc = R0_23_multacc^6./(R0_23_multacc^6 + f*r_23^6);

% the final FRET efficiency is then
% E13 (direct) + E12*E23 (via intermediate)
E = E_13_multacc + E_12_multacc*E_23_multacc;

end