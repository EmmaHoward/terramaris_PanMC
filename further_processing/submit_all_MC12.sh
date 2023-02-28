#!/bin/bash

#
# create postprocessed_outputs data for MC12 
#

year=2007

queue="-p test"
queue="-p short-serial"
#queue="-p short-serial-4hr -A short4hr"

logpath="/home/users/emmah/log/MC12_further/"
scriptpath="/home/users/emmah/python/terramaris_PanMC/further_processing/"

pp="/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/lib/python3.8/site-packages:/home/users/emmah/python/terramaris_PanMC/share:/home/users/emmah/python/emma:/home/users/emmah/python/emma:/home/users/emmah/.conda/envs/stratify/lib/python3.8/site-packages/"

# apparent heat source and moisture sink
sbatch --export=PYTHONPATH=$pp  $queue --array=[1-91] -o $logpath/Q1Q2_${year}_15-%a.o -e $logpath/Q1Q2_${year}-%a.e -t 04:00:00 --mem=24000 $scriptpath/calc_Q1Q2_MC12.py $year

# monthly mean precip diurnal cycle
sbatch $queue  -o $logpath/monprecip_${year}.o -e $logpath/monprecip_${year}.e -t 01:00:00 --mem=24000 extract_diurnal_cycle.py $year MC12

# daily precip
sbatch $queue  -o $logpath/dayprecip_${year}.o -e $logpath/dayprecip_${year}.e -t 01:00:00 --mem=24000 daily_precip.py $year MC12

# Daily mean wind - 850, 200 hPa
sbatch $queue  -o $logpath/wind_${year}.o -e $logpath/wind_${year}.e -t 4:00:00 --mem=38000 regrid_winds_era5.py $year MC12

# sst means
sbatch $queue  --array=[1-3] -o $logpath/sst_%a_${year}.o -e $logpath/sst_%a_${year}.e -t 01:00:00 --mem=24000 calc_sst_means.py $year MC12

# pressure level monthly means
sbatch $queue  -o $logpath/monmean_${year}.o -e $logpath/monmean_${year}.e -t 04:00:00 --mem=24000 monmean_pl_MC12.py $year MC12

