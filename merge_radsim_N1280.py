#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python

import iris
import datetime as dt
from subprocess import check_call
import sys


t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1

scratch = "/work/scratch-pw/emmah/radsim/"
finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/production_runs/%04d%02d_u-cf309/"%(year,(year+1)%100)


def merge(t1,n):
  files = []
  for i in range(1,n*24+1):
     date = t1+dt.timedelta(i/24)
     files.append(scratch+"tma_radsim-himawari_8_ahi-%04d%02d%02d_%02d.nc"%(date.year,date.month,date.day,date.hour))
  data=iris.load(files)
  #data[0].coord("time").units='hours since %04d-11-01 00:00:00'%year
  iris.save(data,scratch+"tma_N1280_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d.nc"%(t1.year,t1.month,t1.day),zlib=True)
  check_call("cp %s/tma_N1280_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d.nc %s/radsim/"%(scratch,t1.year,t1.month,t1.day,finalpath),shell=True)
  for f in files:
    check_call("rm "+f,shell=True)
  check_call("rm %s/tma_N1280_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d.nc"%(scratch,t1.year,t1.month,t1.day),shell=True)

merge(t1,6)
