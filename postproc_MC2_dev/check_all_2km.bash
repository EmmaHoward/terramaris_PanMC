#/bin/bash
set -x


queue="-p short-serial-4hr -A short4hr"
#queue="-p short-serial"

for date in "20151201" "20151207" "20151213" "20151219" "20151225" "20151231" #"20150106" "20150112" "20150118" "20150124" "20150130" "20150205" "20150211" "20150217" "20150223" 
do
  JID_check=`sbatch $queue -o /home/users/emmah/log/MC2_$date.o -e /home/users/emmah/log/MC2_$date.e -t 01:00:00 /home/users/emmah/python/postproc_coupled/check_2km.py $date | cut -d " " -f 4`
done

