#!/bin/bash

year=2007

queue="-p test"
#queue="-p short-serial"
#queue="-p short-serial-4hr -A short4hr"

logpath="/home/users/emmah/log/MC12_further/"
scriptpath="/home/users/emmah/python/terramaris_PanMC/further_processing/"

# apparent heat source and moisture sink
#sbatch $queue --array=[1-91] -o $logpath/Q1Q2_${year}_15-%a.o -e $logpath/Q1Q2_${year}-%a.e -t 01:00:00 --mem=24000 $scriptpath/calc_Q1Q2_MC12.py $year

# monthly mean precip diurnal cycle
sbatch $queue  -o $logpath/monprecip_${year}.o -e $logpath/monprecip_${year}.e -t 01:00:00 --mem=24000 extract_diurnal_cycle.py $year MC12

# sst means
sbatch $queue  --array=[1-3] -o $logpath/monprecip_${year}.o -e $logpath/monprecip_${year}.e -t 01:00:00 --mem=24000 calc_sst_means.py $year MC12

