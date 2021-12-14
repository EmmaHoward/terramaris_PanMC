#/bin/bash
set -x


#queue="-p short-serial-4hr -A short4hr"
queue="-p short-serial"

#for date in "20181201" "20181207" "20181213" "20181219" "20181225" "20181231" "20190106" "20190112" "20190118" "20190124" "20190217" # "20190130" "20190205" "20190211" "20190217" "20190223" 

for date in  "20060130"
do
  JID_check=`sbatch $queue -o /home/users/emmah/log/N1280_$date.o -e /home/users/emmah/log/N1280_$date.e -t 01:00:00 /home/users/emmah/python/terramaris_PanMC/postproc_MC12/check_N1280.py $date | cut -d " " -f 4`
done


