# 2025_ACSNanosciAu_EXAFS
## Related paper
- Journal : ACS Nanoscience Au
- DOI : [10.1021/acsnanoscienceau.5c00013](https://doi.org/10.1021/acsnanoscienceau.5c00013)
- Title : Unravelling Element-Selective Local Structures in Multi-Element Alloy Nanoparticles with EXAFS

## Overview
This repository contains python scripts developed for EXAFS curve fitting analysis of multi-element alloy (high-entropy alloy / compositionally complex alloy) nanoparticles and their application to four samples.
Simultaneous fitting of multiple data at different absorption edges and temperatures is performed while imposing constraints based on physically reasonable and interpretable assumptions.
The details of the procedure is discussed in the article above. The code requires xraylarch version 0.9.65.

## Contents
- scripts : Python scripts developed for the EXAFS analyses
  - tools.py : A collection of codes related to data processing and export
  - run_fit.py : A main script to perform fitting
<br>

- examples : The examples of the structural analyses, The results shown in the above mentioned paper
  - feffinp : The FEFF calculation results (Outputs were not uploaded to reduce contents)
  - in : The experimental data (Athena project files)
  - out : The EXAFS analyses results
    - {samplename}
      - dat : The EXAFS spectra extracted after background subtraction and Fourier transforms (Outputs were not uploaded to reduce contents except for sample "PGM")
      - fit : The EXAFS fitting result of the last cycle of the fitting
      - trial{n} : The results of the {n}th cycle of the fitting
      - samplename_fit.txt : Summary of the fitting result
      - {samplename}\_{n,r,ss}\_{val,err}.csv : The estimated parameters of coordination numbers, atomic radiis, and MSDs
      - {samplename}\_ss\_fit.{txt,png,csv} : The results of Einstein model fitting of MSDs
      - input_pars_{samplename}.csv : A copy of the input parameter file in the same name
      - run_fit.py : A copy of the python script in the same name
  - input_pars_{samplename}.csv : Input files to specify the EXAFS analyses conditions
  - prep_feff.py : A script to perform FEFF calculation before fitting
  - {tools,run_fit}.py : The python scripts
  - seq.ps1 : A script to run a series of fittings

## Usage
Please cite the article above when you publish results using this code.

1. Preparation of input folders/files
  - ./in/{samplename}.prj
    - Athena project file containing experimental data
    - Data in the project should be named as {samplename}\_{temperature}\_{element}
      <br>
      
  - ./feffinp/{element1}{edge}-{element2}/feff0001.dat
    - FEFF calculation results
    - In this version of code, only the first nearest neighbors are considered.
    - If your target has an fcc-based structure, "prep_feff.py" in the "example" is useful.
      <br>
    
  - input_pars_{samplename}.csv
    - Parameters used in the EXAFS analysis
    - List of adjustable parameters
      <br>

      | Parameter | Description |
      | --- | --- |
      | elem | element |
      | edge | absorption edge (K,L3 etc.) |
      | composition | molar fraction of the element |
      | amp | $S_0^2$ parameter used for EXAFS fitting |
      | enot_init | Initial value of $\Delta E_0$ |
      | r_init | Initial value of atomic radius |
      | pre_min | Range used to determine a pre-edge line (min) |
      | pre_max | Range used to determine a pre-edge line (max) |
      | post_min | Range used to determine a post-edge line (min) |
      | post_max | Range used to determine a post-edge line (min) |
      | nnorm | Polynomial order of a post-edge line |
      | bkg_max | Range used to determine a background (max) |
      | rbkg | $R_\mathrm{bkg}$ value used in background subtraction by AUTOBK algorithm |
      | kweight | Exponent of $k$-weighting of EXAFS before forward FT |
      | kmin_ftf | $k$-range of forward FT (min) |
      | kmax_ftf | $k$-range of forward FT (max) |
      | dk_ftf | Tapering parameter of forward FT window (Hanning) |
      | rmin_ftr | $r$-range of reverse FT (min) |
      | rmax_ftr | $r$-range of reverse FT (max) |
      | dr_ftr | Tapering parameter of reverse FT window (Hanning) |

<br>

2. Set parameters in run_fit.py
  - Edit the parameters between "input start" and "input end" comments
<br>

3. Run run_fit.py

## Contact
If you have any question on the scripts, please feel free to contact the authors.
