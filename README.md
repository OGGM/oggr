#Glacier mass balance model of Marzeion et al. 2012 adapted for OGGR 

#Versions

- v1.0 (master_fortran): benchmark version replicating Marzeion 2012 results

##Current (main) branches:

###master:

- using scons to build fortran 
    - to install: e.g. apt-get install scons
    - build main : scons
    - build precompute: scons -Q mmbm_precomp.exe
    - clean main: scons -c
    - clean precompute: scons -c -Q mmbm_precomp.exe

- developed and tested with gfortran

- uses netcdf

- main loop is parallelized with openMP. To use this: 
    - set environment variable: set -x OMP_NUM_THREADS 4
    - set it to your number of CPU-cores as a first estimate
    - enable the two -fopenmp flags in SConstruct
    - BEWARE: if you use parallel version and add new variables, make sure you set the context 
  i.e. private or shared (see beginning of loop)
 
- all python packages are installable via pip or easy_install

- directory structure:
    
    Directories should be symlinked if not wanted in local dir; they are NOT included in the repo

    - _data_ directory
      - {METTYPE}\_anom_interp\_{YEAR}.nc interpolated anomalies
      - {METTYPE}\_t2m_tp_CLIM\_{YEARSTART}_{YEAREND}.nc climatologies
      - cross_validation.nc guess what
    - _precomputed_ directory
      - {REGION}/{METTYPE}\_precomp\_{REGION}_{YEAR}.nc
    - _results_ directory

- NETCDF Conventions:
    - ALL NETCDF files should contain variables in SI units!
    - all lat lons are named latitude, longitude
    - lons range from -180 to + 180
    - time is expressed as decimal years
    - time variable is called month in case of climatology, else time
    - temperature is t2m
    - precip is tp

