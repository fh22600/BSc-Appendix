This directory contains all raw data and Jupyter analysis notebooks for nonequilibrium molecular dynamics simulations of Kemp eliminase variants 1A53-2 (1A0_Sim) and 1A53-2.5 (1A5_Sim).
'miniforge3' was used for the installation of the python package manager miniconda, used to create the python environments to use Jupyter notebooks and execute python scripts. 


Within each directory:
Directories 'Repeat 1-10' contain the raw data from nonequilibrium simulations for each replicate. 
setup.sh scripts were used for setting up directories for the nonequilibrium simulations and removing the substrate from the structures, submit.sh scripts were used for performing MD simulations.
The 'analysis' directory contains all Python scripts and Jupyter notebooks used for analysis of simulation data.
Additional README files accompany additional directories within the 1A0_Sim and 1A5_Sim directories for additional information.

Please note that filepaths used in analysis ntoebooks may need changing as extra directories needed to be made in order to commit data and code to GitHub.
If you have any queries, please contact me at fh22600@bristol.ac.uk. Thank you! 