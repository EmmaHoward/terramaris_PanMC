#!/apps/contrib/jaspy/miniconda_envs/jaspy3.7/m3-4.5.11/envs/jaspy3.7-m3-4.5.11-r20181219/bin/python
import os
import iris
from file_split_2km import file_variables_full,file_variables_pepf
import sys
import numpy as np
import datetime as dt
t0 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

insets,atmos,ocean,radsim = 1,0,0,0

if t0.month >6:
  year = t0.year
else:
  year = t0.year - 1

finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/%04d%02d_u-cc339/"%(year,(year+1)%100)

reinit_step_nc = {"pa":4,"pb":24,"pc":24,"pd":24}   # output, in days

reinit_step_pa = {"density":4,"potential_temperature":4,"pressure_rho_grid":4,"pressure_theta_grid":4,\
                  "upward_air_velocity":4,"eastward_wind":4,"northward_wind":4,"specific_humidity":4,\
                  "radar_reflectivity":24,"cloud_fraction":24,"qcf":24,"qcl":24,"rain":24,"graupel":24}

def check_atmos_insets(stream):
  files = {}
  files_okay = {}
  for var in file_variables_pepf.keys():
    cell_methods,variables = file_variables_pepf[var]
    print("%s: %s"%(stream,var),end=": ")
    expected_length = 6
    files[var] = []
    for i in range(6):
      t = t0+dt.timedelta(i)
      datestring = "%04d%02d%02d"%(t.year,t.month,t.day)
      files[var] += [f for f in os.listdir(finalpath+stream) if (var in f and datestring in f)]
    if len(files[var])!=expected_length:
      print("wrong number of files for %s "%var)
      print("expected: %d, actual: %d"%(expected_length,len(files[var])))
      print("found",files[var])   
      files_okay[var]=False
    else:
      print("correct number of files")
      files_okay[var]=True
  vars_okay = {} 
  for var in file_variables_pepf.keys():
    cell_methods,variables = file_variables_pepf[var]
    reinit = 24
    if files_okay[var]:
      print("%s: %s"%(stream,var),end=": ")
      vars_okay[var]=True
      for f in files[var]:
         data=iris.load(finalpath+stream+"/"+f,variables)
         if len(data) != len(variables):
             print("wrong number of variables for %s:%d"%(f,len(data)))
             vars_okay[var]=False
         for cube in data:
           if len(cube.coords())==0:
              print("cube coord issue on %s at %04d%02d%02d"%(stream,t0.year,t0.month,t0.day))
              print(cube.summary(shorten=True))
              vars_okay[var]=False
           elif (  cube.coord("time").shape[0] != reinit*12): 
              print("cube coord issue at %s"%f)
              print(cube.summary(shorten=True))
              vars_okay[var]=False
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
              print("data corrupted on %s"%f) 
              print(cube.summary(shorten=True))
              vars_okay[var]=False
      if vars_okay[var]:
        print("correct number of variables and data okay")
  print("")
  return (files_okay,vars_okay)


def check_atmos_main():
  files = {}
  files_okay = {}
  for var in file_variables_full.keys():
    stream,cell_methods,variables = file_variables_full[var]
    print("%s: %s"%(stream,var),end=": ")
    reinit = reinit_step_nc[stream]
    if stream=="pa":
      reinit = reinit_step_pa[var]
    expected_length = int(24/reinit)*6
    files[var] = []
    for i in range(6):
      t = t0+dt.timedelta(i)
      datestring = "%04d%02d%02d"%(t.year,t.month,t.day)
      files[var] += [f for f in os.listdir(finalpath+stream) if (var in f and datestring in f)]
    if len(files[var])!=expected_length:
      print("wrong number of files for %s "%var)
      print("expected: %d, actual: %d"%(expected_length,len(files[var])))
      print("found",files[var])   
      files_okay[var]=False
    else:
      print("correct number of files")
      files_okay[var]=True
  vars_okay = {} 
  for var in file_variables_full.keys():
    stream,cell_methods,variables = file_variables_full[var]
    reinit = reinit_step_nc[stream]
    if stream=="pa":
      reinit = reinit_step_pa[var]
    if files_okay[var]:
      print("%s: %s"%(stream,var),end=": ")
      vars_okay[var]=True
      for f in files[var]:
         data=iris.load(finalpath+stream+"/"+f,variables)
         if var == "range":
           if len(data) != 5:
             print("wrong number of variables for %s:%d"%(f,len(data)))
             vars_okay[var]=False
         elif len(data) != len(variables):
             print("wrong number of variables for %s:%d"%(f,len(data)))
             vars_okay[var]=False
         for cube in data:
           if len(cube.coords())==0:
              print("cube coord issue on %s at %04d%02d%02d"%(stream,t0.year,t0.month,t0.day))
              print(cube.summary(shorten=True))
              vars_okay[var]=False
           elif (  (var == "rainfall" and cube.coord("time").shape[0] != reinit*4) \
                 or (var == "daily_mean" and cube.coord("time").shape[0] != 1)\
                 or (var == "integral_agg" and cube.coord("time").shape[0] not in [reinit,8])\
                 or (var in ["surf_agg","range"] and cube.coord("time").shape[0] not in [reinit,1])\
                 or (var not in ["rainfall","daily_mean","integral_agg","surf_agg","range"] and cube.coord("time").shape[0] != reinit)):
              print("cube coord issue at %s"%f)
              print(cube.summary(shorten=True))
              vars_okay[var]=False
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
              print("data corrupted on %s"%f) 
              print(cube.summary(shorten=True))
              vars_okay[var]=False
      if vars_okay[var]:
        print("correct number of variables and data okay")
  print("")
  return (files_okay,vars_okay)


def check_ocean():
  err = []
  t1 = (t0-dt.datetime(year,1,1)).days+1
  dates = [t0 + dt.timedelta(i) for i in range(0,6,1)]
  print("checking kpp for %04d%02d%02d"%(t0.year,t0.month,t0.day))
#  raw = iris.cube.CubeList([cube for cube in raw if cube.name() not in ["surface_altitude","land_binary_mask"]])
  for date in dates:
    if t0.day==1 and t0.month==12:
      rawpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/%04d/%04d%02d%02dT0000Z/"%(year,date.year,date.month,date.day)
    else:
      rawpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/%04d/%04d%02d%02dT0000Z/"%(year,t0.year,t0.month,t0.day)
    raw  = iris.load("%s/KPPocean_%05d.nc"%(rawpath,(date-dt.datetime(year,1,1)).days+1  ))
    processed = iris.load("%s/%s/*%04d%02d%02d.nc"%(finalpath,"kpp",date.year,date.month,date.day))
    processed = iris.cube.CubeList([cube for cube in processed if cube.name() not in ["surface_altitude","land_binary_mask"]])
    if len(raw) != len(processed):
      err.append(1)
      print("mismatch on KPP at %04d%02d%02d: raw contains %d fields, processed contains %d fields"%(t0.year,t0.month,t0.day,len(raw),len(processed)))
  if len(err)==0:
    print("all good!")
  return err



def check_radsim():
  err = []
  for i in range(24):
    t1 = t0+dt.timedelta(i/4)
    print("checking radsim for %04d%02d%02d_%02d"%(t1.year,t1.month,t1.day,t1.hour))
    if not os.path.exists(finalpath+"radsim/tma_2km_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d_%02d.nc"%(t1.year,t1.month,t1.day,t1.hour)):
       print("no radsim file on GWS for %04d%02d%02d"%(t0.year,t0.month,t0.day))
       return [1]
    data=iris.load(finalpath+"radsim/tma_2km_KPPcoupled_bt_himawari_8_ahi-%04d%02d%02d_%02d.nc"%(t1.year,t1.month,t1.day,t1.hour))
    if data[0].shape != (18, 10, 1500, 3200):
       print("radsim issue at %04d%02d%02d"%(t1.year,t1.month,t1.day))
       print(data[0].summary(shorten=True))
       err.append(1)
       try:
         x = data[0][0,0,0,0]
       except RuntimeError:
         print("radsim read issue at %04d%02d%02d"%(t1.year,t1.month,t1.day))
         print(data[0].summary(shorten=True))
         err.append(1)
  if len(err)==0:
    print("radsim okay")
  return err



if atmos:
  atmos_files_okay, atmos_vars_okay = check_atmos_main()
  print("")

if insets:
  pe_files_okay, pe_vars_okay = check_atmos_insets("pe")
  pf_files_okay, pf_vars_okay = check_atmos_insets("pf")
  print("")

if ocean:
  ocean_okay = check_ocean()
  print("")

if radsim:
  radsim_okay = check_radsim()
  print("")

if insets:
  if np.all([pe_files_okay[var] for var in file_variables_pepf]):
    if np.all([pe_vars_okay[var] for var in file_variables_pepf]):
      print("pe okay")
    else:
      print("some pe data issues: see above")
  else:
    print("some pe files missings: see above")
  if np.all([pf_files_okay[var] for var in file_variables_pepf]):
    if np.all([pf_vars_okay[var] for var in file_variables_pepf]):
      print("pf okay")
    else:
      print("some pf data issues: see above")
  else:
    print("some pf files missings: see above")

if atmos:
  if np.all([atmos_files_okay[var] for var in file_variables_full]):
    if np.all([atmos_vars_okay[var] for var in file_variables_full]):
      print("atmos okay")
    else:
      print("some atmos data issues: see above")
  else:
    print("some atmos files missings: see above")

if ocean:
  if len(ocean_okay)==0:
    print("ocean okay")
  else:
    print("some ocean errors: see above")

if radsim:
  if len(radsim_okay)==0:
    print("radsim okay")
  else:
    print("some radsim errors: see above")
