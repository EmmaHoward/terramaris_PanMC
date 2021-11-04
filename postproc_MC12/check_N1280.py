#!/apps/contrib/jaspy/miniconda_envs/jaspy3.7/m3-4.5.11/envs/jaspy3.7-m3-4.5.11-r20181219/bin/python
import os
import datetime as dt
import iris
from subprocess import check_call
import sys
t0 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

if t0.month >6:
  year = t0.year
else:
  year = t0.year - 1

rawpath = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/u-cf309/%04d/"%year
finalpath = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/%04d%02d_u-cf309/"%(year,(year+1)%100)

pp_read_hour = {"pa":0,"pb":0,"pc":0,"pd":12} # input, in hours
reinit_step_nc = {"pa":1,"pb":6,"pc":1,"pd":6}   # output, in days
reinit_step_pp = {"pa":4,"pb":24,"pc":6,"pd":12} # input, in hours

print(t0)

good_to_go = []
err = []


def check(stream,t0):
  err = []
  dates = [t0 + dt.timedelta(i) for i in range(0,6,reinit_step_nc[stream])]
  print("checking %s for %04d%02d%02d"%(stream,t0.year,t0.month,t0.day))
#  raw  = iris.load("%s/%04d%02d%02dT0000Z/tma_%s%04d%02d%02d_%02d.??"%(rawpath,t0.year,t0.month,t0.day,stream,t0.year,t0.month,t0.day,pp_read_hour[stream]))
#  raw = iris.cube.CubeList([cube for cube in raw if cube.name() not in ["surface_altitude","land_binary_mask"]])
  for date in dates:
    try: 
      processed = iris.load("%s/%s/*%04d%02d%02d.nc"%(finalpath,stream,date.year,date.month,date.day))
      processed = iris.cube.CubeList([cube for cube in processed if cube.name() not in ["surface_altitude","land_binary_mask"]])
#      if len(raw) != len(processed):
#        err.append(1)
#        print("mismatch on %s at %04d%02d%02d: raw contains %d fields, processed contains %d fields"%(stream,t0.year,t0.month,t0.day,len(raw),len(processed)))
      for cube in processed:
        if len(cube.coords())==0:
          err.append(1)
          print("cube coord issue on %s at %04d%02d%02d"%(stream,t0.year,t0.month,t0.day))
          print(cube.summary(shorten=True))
        elif cube.coord("time").shape[0] not in [24*reinit_step_nc[stream],4*24*reinit_step_nc[stream],8*reinit_step_nc[stream],reinit_step_nc[stream]]:
          err.append(1)
          print("cube coord issue on %s at %04d%02d%02d"%(stream,t0.year,t0.month,t0.day))
          print(cube.summary(shorten=True))
        try:
          if cube.ndim==4:
            a=cube[0,0,0,0].data
          elif cube.ndim==3:
            a=cube[0,0,0].data
          elif cube.ndim==2:
            a=cube[0,0].data
          else:
            a=cube.data
        except RuntimeError:
          print("data corrupted on %s at %04d%02d%02d"%(stream,t0.year,t0.month,t0.day))
          print(cube.summary(shorten=True))
          err.append(1)
    except OSError:
       err.append(1)
       print("%s at %04d%02d%02d not found"%(stream,t0.year,t0.month,t0.day))
  if len(err)==0:
    print("all good!")
    good_to_go.append("%s_%04d%02d%02d"%(stream,t0.year,t0.month,t0.day))
  return err

def check_ocean(t0):
  err = []
  t1 = (t0-dt.datetime(year,1,1)).days+1
  dates = [t0 + dt.timedelta(i) for i in range(0,6,1)]
  print("checking kpp for %04d%02d%02d"%(t0.year,t0.month,t0.day))
#  raw  = iris.load("%s/%04d%02d%02dT0000Z/KPPocean_%05d.nc"%(rawpath,t0.year,t0.month,t0.day,t1))
#  raw = iris.cube.CubeList([cube for cube in raw if cube.name() not in ["surface_altitude","land_binary_mask"]])
  for date in dates:
    processed = iris.load("%s/%s/*%04d%02d%02d.nc"%(finalpath,"kpp",date.year,date.month,date.day))
    processed = iris.cube.CubeList([cube for cube in processed if cube.name() not in ["surface_altitude","land_binary_mask"]])
#    if len(raw) != len(processed):
#      err.append(1)
#      print("mismatch on KPP at %04d%02d%02d: raw contains %d fields, processed contains %d fields"%(t0.year,t0.month,t0.day,len(raw),len(processed)))
  if len(err)==0:
    print("all good!")
    good_to_go.append("%s_%04d%02d%02d"%("kpp",t0.year,t0.month,t0.day))
  return err


def check_radsim(t0):
  err = []
  for i in range(6):
    t1 = t0+dt.timedelta(i)
    if not os.path.exists(finalpath+"radsim/tma_N1280_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d.nc"%(t1.year,t1.month,t1.day)):
       print("no radsim file on GWS for %04d%02d%02d"%(t0.year,t0.month,t0.day))
       return [1]
    data=iris.load(finalpath+"radsim/tma_N1280_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d.nc"%(t1.year,t1.month,t1.day))
    print("checking radsim for %04d%02d%02d"%(t1.year,t1.month,t1.day))
    if data[0].shape != (72, 10, 428, 538):
       print("radsim issue at %04d%02d%02d"%(t1.year,t1.month,t1.day))
       print(data[0].summary(shorten=True))
       err.append(1)
    else:
       print("all good!")
       good_to_go.append("%s_%04d%02d%02d"%("radsim",t1.year,t1.month,t1.day))
  return err

def remove(stream,t0):
  dates = [t0 + dt.timedelta(i/24) for i in range(0,6*24,reinit_step_pp[stream])]
  for date in dates:
    cmd  = "ls %s/%04d%02d%02dT0000Z/tma_%s%04d%02d%02d_%02d.pp"%(rawpath,t0.year,t0.month,t0.day,stream,date.year,date.month,date.day,date.hour)
    #print(cmd)
    check_call(cmd,shell=True)

for stream in ["pa","pb","pc","pd"][:]:
  err += check(stream,t0)
err += check_radsim(t0)
err += check_ocean(t0)

print(err)
print(good_to_go)


