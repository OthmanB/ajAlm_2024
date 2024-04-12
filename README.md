# What is this
This is the ensemble of programs required to analyse MCMC results related to the
project of analysing all of the Kepler stars in a systematic manner to:
   - identify which model complexity is the statistical significant (based on the global likelihood computation)
   - prepare plots for all aspects of rotational parameters (a1,a2,a3,a4) and/or activity-related parameters (epsilon, theta0, delta). See Benomar et al. 2023 for definitions
   - establish graphs and plots that show active latitudes for the ensemble of stars analysed in this large scale analysis


# Instructions and dependencies
This is still work in progress. Instructions will be given once the code is consolidated.
However, most of the code here is in python3. Some scripts are in bash as well, but these are supposed to be used in the HPC and/or on the raw MCMC results
python programs often need to call binaries of the TAMCMC program. The version used here is the 1.86.78, which is available as Linux-based x86 binaries at https://github.com/OthmanB/TAMCMC-C/releases/tag/v.186.78 . For MacOS, please compile using the source code and cmake. The program is untested in Windows.

Main data are currently stored in my NAS. These weight a few Tb. So obviously they cannot get in github. 
Eventually, lightweight statistical products will be made public. And main data may be made available through a SQL database.


# Directory structure
   - programs: All of the programs
   - External_data: Ensemble of the papers and data sources used to create tables of non-seismic parameters

