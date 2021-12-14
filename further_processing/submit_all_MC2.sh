#!/bin/bash

year=2016

queue="-p test"
#queue="-p short-serial"
#queue="-p short-serial-4hr -A short4hr"

logpath="/home/users/emmah/log/MC2_further/"
scriptpath="/home/users/emmah/python/terramaris_PanMC/further_processing/"

# apparent heat source and moisture sink
#sbatch $queue --array=[1-91] -o $logpath/rho_${year}-%a.o -e $logpath/rho_${year}-%a.e -t 04:00:00 --mem=128000 $scriptpath/coarsen_rho_MC2.py $year
#sbatch $queue --array=[1-91] -o $logpath/Q1Q2_${year}-%a.o -e $logpath/Q1Q2_${year}-%a.e -t 04:00:00 --mem=48000 $scriptpath/calc_Q1Q2_MC2.py $year

# monthly mean precip diurnal cycle
#sbatch $queue  -o $logpath/monprecip_${year}.o -e $logpath/monprecip_${year}.e -t 01:00:00 --mem=24000 extract_diurnal_cycle.py $year MC2

# daily precip on MC12 grid
#sbatch $queue --array=[1-3] -o $logpath/dayprecip_${year}_%a.o -e $logpath/dayprecip_${year}_%a.e -t 01:00:00 --mem=24000 coarsen_precip_MC2.py $year

# sst means
#sbatch $queue --array=[1-3] -o $logpath/monsst_${year}_%a.o -e $logpath/monsst_${year}_%a.e -t 01:00:00 --mem=24000 calc_sst_means.py $year MC2

sbatch $queue  -o $logpath/wind_${year}.o -e $logpath/wind_${year}.e -t 04:00:00 --mem=128000 regrid_winds_era5.py $year MC2
