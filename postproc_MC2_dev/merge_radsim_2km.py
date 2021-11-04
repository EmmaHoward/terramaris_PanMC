#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python

import iris
import datetime as dt
from subprocess import check_call
import sys
import os


t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1

scratch = "/work/scratch-pw/emmah/radsim_2km/"
finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/%04d%02d_u-cc339/"%(year,(year+1)%100)


def merge(t1,n):
  allfiles = []
  for i in range(n*24+1):
     date = t1+dt.timedelta(i/24)
     for minute in [0,20,40]:
       allfiles.append(scratch+"tm2a_radsim-himawari_8_ahi-%04d%02d%02d_%02d%02d.nc"%(date.year,date.month,date.day,date.hour,minute))
  allfiles = allfiles[1:-2]
  for i in range(4):
    files = allfiles[i*18:(i+1)*18]
    hour = i*6
    data=iris.load(files)
    print(data)
    #data[0].coord("time").units='hours since %04d-11-01 00:00:00'%year
    iris.save(data,scratch+"tma_2km_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d_%02d.nc"%(t1.year,t1.month,t1.day,hour),zlib=True)
    check_call("cp %s/tma_2km_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d_%02d.nc %s/radsim/"%(scratch,t1.year,t1.month,t1.day,hour,finalpath),shell=True)
    for f in files:
      check_call("rm "+f,shell=True)
    check_call("rm %s/tma_2km_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d_%02d.nc"%(scratch,t1.year,t1.month,t1.day,hour),shell=True)

job = int(os.environ["SLURM_ARRAY_TASK_ID"])-1
merge(t1+dt.timedelta(job),1)
