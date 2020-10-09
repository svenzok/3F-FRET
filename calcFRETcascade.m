function E = calcFRETcascade(r_12,r_13,r_23,R0_12,R0_13,R0_23)
%CALCFRETCASCADE calculates cascading FRET efficiencies
%
%   inputs:
%   |r_12|, distance between F1 and F2 (nm)
%   |r_23|, distance between F2 and F3 (nm)
%   |r_13|, distance between F1 and F3 (nm)
%   |R0_12|, Förster distance between F1 and F2 (nm)
%   |R0_23|, Förster distance between F2 and F3 (nm)
%   |R0_13|, Förster distance between F1 and F3 (nm)
%
%   output:
%   |E|, energy transfer efficiency from F1 to F3, by summing up the
%   possible pathways of either 1->3 (direct) or 1->2->3 (via the
%   intermediate F2)

% The following calculations are based on:
% Sun et al., Three-Color Spectral FRET Microscopy Localizes Three
% Interacting Proteins in Living Cells,
% Biophys. J. 2010, 99, 4, 1274-1283
% https://doi.org/10.1016/j.bpj.2010.06.004
E_12 = (R0_12.*r_13)^6/((R0_12.*r_13)^6+(R0_13.*r_12)^6+(r_12.*r_13)^6);
E_13 = (R0_13.*r_12)^6/((R0_12.*r_13)^6+(R0_13.*r_12)^6+(r_12.*r_13)^6);
E_23 = R0_23^6/(R0_23^6 + r_23^6);

% the final FRET efficiency is then E_13 (direct) + E_12 * E_23 (via the
% intermediate)
E = E_13 + E_12 * E_23;

end