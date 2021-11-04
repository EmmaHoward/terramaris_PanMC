#!/bin/bash
set -x 
  
date=$1


queue="-p short-serial-4hr -A short4hr"
#queue="-p short-serial"

  

#sbatch $queue -J MC2_pa1 --array=[1,6-10] --mem=96000 -o /home/users/emmah/log/postproc_2km/pp1_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp1_%a_${date}.e -t 04:00:00 postproc_2km.py $date 0 6
#sbatch $queue -J MC2_pa2 --array=[2-5,11-14] --mem=24000 -o /home/users/emmah/log/postproc_2km/pp1_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp1_%a_${date}.e -t 04:00:00 postproc_2km.py $date 0 3
#sbatch $queue -J MC2_pa3 --array=[2-5,11-14] --mem=24000 -o /home/users/emmah/log/postproc_2km/pp1_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp1_%a_${date}.e -t 04:00:00 postproc_2km.py $date 3 6
#sbatch $queue -J MC2_pb  --array=[15-21] --mem=24000 -o /home/users/emmah/log/postproc_2km/pp1_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp1_%a_${date}.e  -t 03:00:00 postproc_2km.py $date 0 6
#for i in $(seq 0 5)
#  do
#  sbatch $queue -J MC2_pc${i} --array=[22-32] --mem=96000 -o /home/users/emmah/log/postproc_2km/pp1_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp1_%a_${date}.e  -t 04:00:00 postproc_2km.py $date $i $((i+1))
#  sbatch $queue --array=[2] -o /home/users/emmah/log/postproc_2km/pp2_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp2_%a_${date}.e -t 04:00:00 --mem=64000 postproc_2km_insets_ff.py $date $i $((i+1)) pe
#  sbatch $queue --array=[2] -o /home/users/emmah/log/postproc_2km/pp2_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp2_%a_${date}.e -t 04:00:00 --mem=64000 postproc_2km_insets_ff.py $date $i $((i+1)) pf
#  done
#sbatch $queue -J MC2_pd --array=[33-43] --mem=30000 -o /home/users/emmah/log/postproc_2km/pp1_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp1_%a_${date}.e  -t 04:00:00 postproc_2km.py $date 0 6

#sbatch $queue -J MC2_kpp --array=[1-11] --job-name="pp_kpp_2km" -o /home/users/emmah/log/postproc_2km/kpp_${date}_%a.o -e /home/users/emmah/log/postproc_2km/kpp_${date}_%a.e  -t 04:00:00 --mem=24000 /home/users/emmah/python/postproc_coupled/postproc_2km_ocean.py ${date}

#sbatch $queue --array=[12] -o /home/users/emmah/log/postproc_2km/pp2_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp2_%a_${date}.e -t 04:00:00 --mem=64000 postproc_2km_insets_ff.py $date 0 6 pe
#sbatch $queue --array=[12] -o /home/users/emmah/log/postproc_2km/pp2_%a_${date}.o -e /home/users/emmah/log/postproc_2km/pp2_%a_${date}.e -t 04:00:00 --mem=64000 postproc_2km_insets_ff.py $date 0 6 pf

#sbatch $queue --array=[1-36] -o /home/users/emmah/log/postproc_2km/tar_%a_${date}.o -e /home/users/emmah/log/postproc_2km/tar_%a_${date}.e  -t 01:00:00 untar_pepf $date

#JID_radsim=`sbatch -J MC2_radsim $queue --array=[1-144] -o /home/users/emmah/log/radsim_2km/radsim_%a_${date}.o -e /home/users/emmah/log/radsim_2km/radsim_%a_${date}.e  -t 04:00:00 --mem=24000 run_radsim_2km.py ${date}  | cut -d " " -f 4`

sbatch $queue --array=[1-4] -o /home/users/emmah/log/radsim_2km/radsim_${date}.o -e /home/users/emmah/log/radsim_2km/radsim_${date}.e -t 01:00:00 /home/users/emmah/python/postproc_coupled/merge_radsim_2km.py ${date}
#sbatch $queue --array=[1-6] -o /home/users/emmah/log/radsim_2km/radsim_${date}.o -e /home/users/emmah/log/radsim_2km/radsim_${date}.e --dependency=afterok:$JID_radsim -t 01:00:00 /home/users/emmah/python/postproc_coupled/merge_radsim_2km.py ${date}

