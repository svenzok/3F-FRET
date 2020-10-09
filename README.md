# 3F-FRET
3F-FRET contains calculation scripts for the publication ["Three-fluorophore Förster resonance energy transfer for the analysis of ternary protein association in plant cells"](https://www.biorxiv.org/content/10.1101/722124v1), and a simple linear unmixing script to separate three single spectral components from a sum spectrum.

### original spectra and photophysical parameters
The spectra and photophysical parameters of mTurquoise2, mVenus and mRFP1 that were used in the above publication are stored in
- FP_original_spectra.mat
- FP_photophysical_params.mat
using the original published data for these fluorescent proteins.

### calculation of Förster distances and 10% FRET distances
The data in Table 1 of the above publication were calculated using
- three_fluorophore_FRET_Table1.m

### spectral unmixing
The main program
- unmix_tripleFRET

applies linear unmixing to a measured spectrum to estimate the ratio of the underlying three components (mTurquoise2, mVenus and mRFP1 for the above publication, but any viable combination can be processed here). This program also produces a figure to visualize the relative contribution of the components to the sum spectrum.

### sub-programs
The following programs are needed for full functionality of the above main scripts, but might also be used in another context (see description in the respective programs): 
- interpol_spec.m
- fit_spec_em.m
- fit_spec_em_logn.m
- align_spec.m
- calc_R0.m
- calcFRETcascade.m
- calc_R0_multacc.m
- calcFRETcascade_multacc.m

# Requirements
- MATLAB R2020a or newer
  - Curve Fitting Toolbox

# Copyright and Software License
Copyright (c) SFB1101, ZMPB and IPTC, University of Tübingen, Tübingen.

The scripts in 3F-FRET are licensed under the [GNU GPL](https://www.gnu.org/licenses/)

# How to cite 3F-FRET
If you use any of the scripts in 3F-FRET to process your data, please, cite our [paper](https://www.biorxiv.org/content/10.1101/722124v1):
- Nina Glöckner, Sven zur Oven-Krockhausa, Leander Rohr, Frank Wackenhut, Moritz Burmeister, Friederike Wanke, Eleonore Holzwart, Alfred J. Meixner, Sebastian Wolf and Klaus Harter. Three-fluorophore Förster resonance energy transfer for the analysis of ternary protein association in plant cells. bioRxiv preprint doi: [https://doi.org/10.1101/722124](https://doi.org/10.1101/722124)
