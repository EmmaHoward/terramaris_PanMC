import iris
import iris.plot as iplt
from panMC import panMC
import cmocean
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import subprocess
from iris.coord_categorisation import add_hour,add_month
from iris.util import equalise_attributes
import precip_bias
from scipy.ndimage import convolve
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/precip/"
MC2_path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"
ref_path = "/gws/nopw/j04/klingaman/datasets/GPM_IMERG/monthly/"
P12 = iris.cube.CubeList()
Pref = iris.cube.CubeList()
cx = iris.Constraint(longitude = lambda x: 84 < x < 161)
cy = iris.Constraint(latitude  =  lambda y: -21 < y < 21)

def plot(MC12years,MC2years,figname,fig):
  P12 = iris.cube.CubeList()
  P2 = iris.cube.CubeList()
  Pref = iris.cube.CubeList()
  for year in MC12years:
    print(year)
    P12.append(precip_bias.load(year,"MC12",MC12_path).collapsed("time",iris.analysis.MEAN))
    Pref.append(precip_bias.GPM(year,P12[-1]).collapsed("time",iris.analysis.MEAN))
    P12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    Pref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    if year in MC2years:
      tmp = iris.load([MC2_path+"%04d%02d_hourly_coarsened_rainfall.nc"%(yr,mon) for (yr,mon) in [(year,12),(year+1,1),(year+1,2)]])
      iris.util.equalise_attributes(tmp)
      for cube in tmp:
         cube.coord('time').convert_units("hours since 2003-11-01") 
      tmp = tmp.concatenate_cube().collapsed("time",iris.analysis.MEAN)*24*4
      P2.append(tmp)
      P2[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  equalise_attributes(P2)
  P2=P2.merge_cube()
  P12=P12.merge_cube()
  P2=P2.collapsed(["longitude"],iris.analysis.MEAN)
  P12=P12.collapsed(["longitude"],iris.analysis.MEAN)
  Pref=Pref.merge_cube()
  Pref=Pref.collapsed(["longitude"],iris.analysis.MEAN)

  P2.data = convolve(P2.data,np.ones((1,15))/15,mode="nearest")
  P12.data = convolve(P12.data,np.ones((1,15))/15,mode="nearest")
  Pref.data = convolve(Pref.data,np.ones((1,15))/15,mode="nearest")

#  P2 = P2.rolling_window("latitude",iris.analysis.MEAN,10)
#  P12 = P12.rolling_window("latitude",iris.analysis.MEAN,10)
#  Pref = Pref.rolling_window("latitude",iris.analysis.MEAN,10) 
  lines = ["-","--",":","-."]
  plt.subplot(321)
  plt.title("MC2")
  if len(MC2years)==1:
      iplt.plot(P2.coord("latitude")[6:-6],P2[6:-6],lw=1)#,label="%04d-%02d"%(year,(year+1)%100))
  else:
    for i,year in enumerate(MC2years):
      iplt.plot(P2.coord("latitude")[6:-6],P2[i,6:-6],lw=[1,1,2,1][i%4],ls = lines[i%4])#,label="%04d-%02d"%(year,(year+1)%100))
  plt.xlim(-20,20)
  plt.plot([-15,-15,15,15,-15],[-1,16,16,-1,-1],"k",lw=1)
  plt.ylabel("Precip (mm/day)")
  plt.ylim(0,15)
  plt.grid()

  plt.subplot(323)
  plt.plot([-15,-15,15,15,-15],[-1,16,16,-1,-1],"k",lw=1)
  plt.title("MC12")
  plt.ylabel("Precip (mm/day)")
  plt.ylim(0,15)
  plt.xlim(-20,20)
  plt.grid()
  for i,year in enumerate(MC12years):
    iplt.plot(P12.coord("latitude")[6:-6],P12[i,6:-6],lw=[1,1,2,1][i%4],ls = lines[i%4])#,label="%04d-%02d"%(year,(year+1)%100))

  plt.subplot(325)
  plt.plot([-15,-15,15,15,-15],[-1,16,16,-1,-1],"k",lw=1)
  plt.xlim(-20,20)
  plt.ylim(0,15)
  plt.title("GPM-IMERG")
  plt.ylabel("Precip (mm/day)")
  plt.xlabel("Latitude")
  plt.grid()
  for i,year in enumerate(MC12years):
    iplt.plot(Pref.coord("latitude"),Pref[i],label="%04d-%02d"%(year,(year+1)%100),lw=[1,1,2,1][i%4],ls = lines[i%4])

  fig.subplots_adjust(top =0.95,bottom=0.18,left=0.08,right=0.98,hspace=0.3,wspace=0.06)
  fig.legend(loc="lower left",ncol=4)
  #fig.savefig(figname)
  #plt.show()
