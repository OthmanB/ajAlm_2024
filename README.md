# What is this
This is the ensemble of programs required to analyse MCMC results related to the
project of analysing all of the Kepler stars in a systematic manner to:
   - identify which model complexity is the statistical significant (based on the global likelihood computation)
   - prepare plots for all aspects of rotational parameters (a1,a2,a3,a4) and/or activity-related parameters (epsilon, theta0, delta). See Benomar et al. 2023 for definitions
   - establish graphs and plots that show active latitudes for the ensemble of stars analysed in this large scale analysis

The model that are handled (and for which MCMC has been performed in NYUAD HPC) are encoded as a sequence of 4 numbers abcd.
The encoding is as follow:
  - a : Relates to a1. Can take values 0 (no rotation) or 1 (using a1)
  - b : Relates to aspherity. Can take values:
       - 0 : no activity accounted for
       - 1 : using a2 and a4 to describe the stellar asphericity, hence the activity + centrifugal effect
       - 2 : using Ajlm(theta0, delta) to describe the activity + centrifugal effect. These contain submodels:
               - Gate (G in short) : Assume the active zone being described by a gate function.
               - Triangle (T in short) : Assume that the active zone is being described by a triangle function. This is shown to very well match the Sun long term activity.
             Each of these submodels are declined into three possible complexity:
               - decompose_-1: The Ajlm is not decomposed into a-coefficients, meaning that if the data contain l=0,1,2,3, all these will be accounted for. 
               - decompose_1 : The Ajlm is decomposed into a-coefficients. Only effects of parity coefficient a2,a4 are accounted for (we drop a6). This allows to ensure consistency among all of the explored models. In particular tests showed that adding more than a1 to l=3 makes the fit unstable. This aleviate that issue
               - decompose_2 : The Ajlm is decomposed into a-coefficients. Only effects of parity coefficient a2 are accounted for (we drop a4 and a6). This allows to account for the simplest activity case.

  - c : Relates to the latitudinal differential rotation. Can take values 0 (no differential rotation) or 1 (using a3)
  - d : Relates to the Lorentzian asymetry. Can take values 0 (no asymetry) or 1 (using model of Benomar et al. 2017)
Models that we computed are:
  - 1001 : a1 with asymetry
  - 1101 : a1,a2,a4 with asymetry
  - 1111 : a1,a2,a4 with asymetry
  - 1201 : a1,Ajlm (all submodels and complexity) with asymetry
  - 1211 : a1,Ajlm (all submodels and complexity), a3 with asymetry

Nearly 1000 combinations are in total avaiable. Totaling a compute time of ~ 2 weeks x 1000 cpus

# Instructions and dependencies
This is still work in progress. Instructions will be given once the code is consolidated.
However, most of the code here is in python3. Some scripts are in bash as well, but these are supposed to be used in the HPC and/or on the raw MCMC results
python programs often need to call binaries of the TAMCMC program. The version used here is the 1.86.78, which is available as Linux-based x86 binaries at https://github.com/OthmanB/TAMCMC-C/releases/tag/v.186.78 . For MacOS, please compile using the source code and cmake. The program is untested in Windows.

Main data are currently stored in my NAS. These weight a few Tb. So obviously they cannot get in github. 
Eventually, lightweight statistical products will be made public. And main data may be made available through a SQL database.


# Directory structure
   - programs: All of the programs
   - External_data: Ensemble of the papers and data sources used to create tables of non-seismic parameters

