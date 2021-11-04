#/bin/bash
set -x

date=20151225

#JID_atmos=`sbatch -p short-serial-4hr -A short4hr --array=[1-45] --job-name="pp_atmos_N1280" -o /home/users/emmah/log/postproc_N1280/pp1_$date_%a.o -e /home/users/emmah/log/postproc_N1280/pp1_$date_%a.e  -t 04:00:00 --mem=24000 /home/users/emmah/python/postproc_coupled/postproc_N1280.py $date | cut -d " " -f 4`

#JID_ocean=`sbatch -p short-serial-4hr -A short4hr --array=[1-11] --job-name="pp_kpp_N1280" -o /home/users/emmah/log/postproc_N1280/kpp_$date_%a.o -e /home/users/emmah/log/postproc_N1280/kpp_$date_%a.e  -t 04:00:00 --mem=24000 /home/users/emmah/python/postproc_coupled/postproc_N1280_ocean.py $date | cut -d " " -f 4`

#JID_radsim=`sbatch -p short-serial-4hr -A short4hr --array=[1-36] -o /home/users/emmah/log/radsim_N1280/radsim_$date_%a.o -e /home/users/emmah/log/radsim_N1280/radsim_$date_%a.e -t 01:00:00 /home/users/emmah/python/postproc_coupled/run_radsim_N1280.py $date | cut -d " " -f 4`  

#JID_merge_radsim=`sbatch -p short-serial-4hr -A short4hr -o /home/users/emmah/log/radsim_N1280/radsim_$date.o -e /home/users/emmah/log/radsim_N1280/radsim_$date.e --dependency=afterok:$JID_radsim -t 01:00:00 /home/users/emmah/python/postproc_coupled/merge_radsim_N1280.py $date | cut -d " " -f 4`  

#JID_check=`ssh sci2 "sbatch -p short-serial-4hr -A short4hr -o /home/users/emmah/log/N1280_$date.o -e /home/users/emmah/log/N1280_$date.e --dependency=afterok:$JID_merge_radsim:$JID_ocean:$JID_atmos -t 01:00:00 /home/users/emmah/python/postproc_coupled/check_N1280.py $date | cut -d ' ' -f 4"` 
JID_check=`ssh sci2 "sbatch -p short-serial-4hr -A short4hr -o /home/users/emmah/log/N1280_$date.o -e /home/users/emmah/log/N1280_$date.e -t 01:00:00 /home/users/emmah/python/postproc_coupled/check_N1280.py $date | cut -d ' ' -f 4"`

echo $JID_check 

