
import os
import iris
from file_split_2km import file_variables_full,file_variables_pepf
import sys
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
t0 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

def nsp(i):
  a=int(np.ceil(np.sqrt(i)))
  b=int(np.ceil(i/a))
  return a,b


if t0.month >6:
  year = t0.year
else:
  year = t0.year - 1

finalpath = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/%04d%02d_u-cc339/"%(year,(year+1)%100)

reinit_step_nc = {"pa":4,"pb":24,"pc":24,"pd":24}   # output, in days

reinit_step_pa = {"density":4,"potential_temperature":4,"pressure_rho_grid":4,"pressure_theta_grid":4,\
                  "upward_air_velocity":4,"eastward_wind":4,"northward_wind":4,"specific_humidity":4,\
                  "radar_reflectivity":24,"cloud_fraction":24,"qcf":24,"qcl":24,"rain":24,"graupel":24}

fig=plt.figure(figsize=(12,12))
def check_atmos_main():
  for j,var in enumerate(file_variables_full.keys()):
    stream,cell_methods,variables = file_variables_full[var]
    print("%s: %s"%(stream,var),end=": ")
    reinit = reinit_step_nc[stream]
    if stream=="pa":
      reinit = reinit_step_pa[var]
    files = []
    for i in range(6):
      t = t0+dt.timedelta(i)
      datestring = "%04d%02d%02d"%(t.year,t.month,t.day)
      files += [f for f in os.listdir(finalpath+stream) if (var in f and datestring in f)]
    f = files[np.random.randint(0,len(files))]
    data=iris.load(finalpath+stream+"/"+f,variables)
    for j1,cube in enumerate(data):
      if stream == "pa":
        ax = plt.subplot(3,5,j+1+j1,projection=ccrs.PlateCarree())
      if stream == "pb":
        ax = plt.subplot(3,5,1+j1%14,projection=ccrs.PlateCarree())
      if stream == "pc":
        ax = plt.subplot(3,4,j-20+j1,projection=ccrs.PlateCarree())
      if stream == "pd":
        ax = plt.subplot(3,5,j-31+j1,projection=ccrs.PlateCarree())
      ax.coastlines()
      if cube.ndim == 4:
          nt,nz,ny,nx = cube.shape
          i1,i2 = np.random.randint(0,nt),np.random.randint(0,nz)
          sub = cube[i1,i2]
          name=" %d %d"%(i1,i2)
      elif cube.ndim == 3:
          nt,ny,nx = cube.shape
          i1 = np.random.randint(0,nt)
          sub = cube[i1]
          name=" %d"%(i1)
      elif cube.ndim ==2:
          sub = cube
          name=" full"
      v= sub.collapsed(["latitude","longitude"],iris.analysis.PERCENTILE,percent=(5,95)).data
      a=iplt.pcolormesh(sub,vmin=v[0],vmax=v[1], axes=ax)
      plt.colorbar(a,ax=ax,orientation="horizontal")
      plt.title(cube.name()+name,fontsize="small")
      if var=="surface_inst" and j1==13:
         plt.savefig("check_%s_%s_1.png"%(stream,var))
         fig.clf()
         fig.set_figwidth(12)
         fig.set_figheight(12)
    if stream =="pb":
       plt.savefig("check_%s_%s.png"%(stream,var))
       fig.clf()
       fig.set_figwidth(12)
       fig.set_figheight(12)
    print(j)
    if j in [13,20,31,42]:
       plt.savefig("check_%s.png"%stream)
       fig.clf()
       fig.set_figwidth(12)
       fig.set_figheight(12)


def check_atmos_inset(stream):
  j1=0
  for var in file_variables_pepf.keys():
    cell_methods,variables = file_variables_pepf[var]
    files = []
    for i in range(6):
      t = t0+dt.timedelta(i)
      datestring = "%04d%02d%02d"%(t.year,t.month,t.day)
      files += [f for f in os.listdir(finalpath+stream) if (var in f and datestring in f)]
    f = files[np.random.randint(0,len(files))]
    data=iris.load(finalpath+stream+"/"+f,variables)
    for j,cube in enumerate(data):
      if len(data)==1:
        ax = plt.subplot(2,4,j+1+j1,projection=ccrs.PlateCarree())
        j1+=1
      elif len(data)>1:
        a,b = nsp(len(data))
        ax = plt.subplot(a,b,j+1,projection=ccrs.PlateCarree())
      ax.coastlines()
      if cube.ndim == 4:
          nt,nz,ny,nx = cube.shape
          i1,i2 = np.random.randint(0,nt),np.random.randint(0,nz)
          sub = cube[i1,i2]
          name=" %d %d"%(i1,i2)
      elif cube.ndim == 3:
          nt,ny,nx = cube.shape
          i1 = np.random.randint(0,nt)
          sub = cube[i1]
          name=" %d"%(i1)
      elif cube.ndim ==2:
          sub = cube
          name=" full"
      a=iplt.pcolormesh(sub,axes=ax)
      plt.colorbar(a,ax=ax,orientation="horizontal")
      plt.title(cube.name()+name,fontsize="small")
    if var == "upward_air_velocity" or len(data)>1:
       plt.savefig("check_%s_%s.png"%(stream,var))
       fig.clf()
       fig.set_figwidth(12)
       fig.set_figheight(12)



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
  print("checking radsim for %04d%02d%02d"%(t0.year,t0.month,t0.day))
  for i in range(24):
    t1 = t0+dt.timedelta(i/4)
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

#check_atmos_main()
check_atmos_inset("pe")
