#!/bin/bash

#
# create postprocessed_outputs data for MC2 
#

set -x
year=2012

#queue="-p test"
#queue="-p short-serial"
queue="-p short-serial-4hr -A short4hr"

logpath="/home/users/emmah/log/MC2_further/"
scriptpath="/home/users/emmah/python/terramaris_PanMC/further_processing/"

pp="/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/lib/python3.8/site-packages:/home/users/emmah/python/terramaris_PanMC/share:/home/users/emmah/python/emma:/home/users/emmah/python/emma:/home/users/emmah/.conda/envs/stratify/lib/python3.8/site-packages/"


# monthly mean precip diurnal cycle
sbatch $queue  --export=PYTHONPATH=$pp    -o $logpath/monprecip_${year}.o -e $logpath/monprecip_${year}.e -t 04:00:00 --mem=24000 extract_diurnal_cycle.py $year MC2

# daily precip on MC12 grid
sbatch $queue  --export=PYTHONPATH=$pp   --array=[1-3] -o $logpath/dayprecip_${year}_%a.o -e $logpath/dayprecip_${year}_%a.e -t 01:00:00 --mem=24000 coarsen_precip_MC2.py $year

# sst means
sbatch $queue  --export=PYTHONPATH=$pp   --array=[1-3] -o $logpath/monsst_${year}_%a.o -e $logpath/monsst_${year}_%a.e -t 01:00:00 --mem=24000 calc_sst_means.py $year MC2

# 850 and 200 hPa daily winds on era-5 grid
sbatch -p short-serial  --export=PYTHONPATH=$pp    -o $logpath/wind_${year}.o -e $logpath/wind_${year}.e -t 10:00:00 --mem=128000 regrid_winds_era5.py $year MC2

# monthly mean pressure levels
JID=`sbatch $queue  --export=PYTHONPATH=$pp   --array=[1-66] -o $logpath/monmean_${year}_%a.o -e $logpath/monmean_${year}_%a.e -t 4:00:00 --mem=128000 monmean_pl_MC2.py $year MC2 calc  | cut -d " " -f 4""`
sbatch --export=PYTHONPATH=$pp  $queue -o $logpath/monmean_${year}.o -e $logpath/monmean_${year}.e -t 4:00:00 --dependency=afterok:$JID  --mem=128000 monmean_pl_MC2.py $year MC2 collate 

# apparent heat source and moisture sink
JID=`sbatch --export=PYTHONPATH=$pp   $queue --array=[1-90] -o $logpath/rho_${year}-%a.o -e $logpath/rho_${year}-%a.e -t 04:00:00 --mem=128000 $scriptpath/coarsen_rho_MC2.py $year  | cut -d " " -f 4""`
sbatch --export=PYTHONPATH=$pp   $queue --array=[1-90] -o $logpath/Q1Q2_${year}-%a.o -e $logpath/Q1Q2_${year}-%a.e --dependency=afterok:$JID -t 04:00:00 --mem=48000 $scriptpath/calc_Q1Q2_MC2.py $year
#sbatch --export=PYTHONPATH=$pp   $queue --array=[1-90] -o $logpath/Q1Q2_${year}-%a.o -e $logpath/Q1Q2_${year}-%a.e  -t 04:00:00 --mem=48000 $scriptpath/calc_Q1Q2_MC2.py $year


