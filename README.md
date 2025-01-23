# 2025_ACSNanosciAu_EXAFS
## Related paper
- Journal : ACS Nanoscience Au
- DOI :
- Title : Unravelling Element-Selective Local Structures in Multi-Element Alloy Nanoparticles with EXAFS

## Contents
- scripts : Python scripts developed for the EXAFS analyses
  - tools.py : A collection of codes related to data processing and export
  - run_fit.py : A main script to perform fitting
<br>

- examples : The examples of the structural analyses, The results shown in the above mentioned paper
  - feffinp : A folder where the FEFF input and outputs are exported
  - in : The experimental data (Athena project files)
  - out : The EXAFS analyses results
    - {samplename}
      - dat : The EXAFS spectra extracted after background subtraction and Fourier transforms
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
  - seq.ps1 : The executive file of a series of fitting
