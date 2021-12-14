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
from iris.experimental.equalise_cubes import equalise_attributes
import precip_bias
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monthly_diurnal_precip_amounts/"
ref_path = "/gws/nopw/j04/klingaman/datasets/GPM_IMERG/monthly/"
P12 = iris.cube.CubeList()
Pref = iris.cube.CubeList()
cx = iris.Constraint(longitude = lambda x: 84 < x < 161)
cy = iris.Constraint(latitude  =  lambda y: -21 < y < 21)

def plot(years,figname):
  P12 = iris.cube.CubeList()
  Pref = iris.cube.CubeList()
  for year in years:
    print(year)
    P12.append(precip_bias.load(year,"MC12",MC12_path).collapsed("time",iris.analysis.MEAN))
    Pref.append(precip_bias.GPM(year,P12[-1]).collapsed("time",iris.analysis.MEAN))
    P12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    Pref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  
  P12=P12.merge_cube()
  P12=P12.collapsed(["longitude"],iris.analysis.MEAN)
  Pref=Pref.merge_cube()
  Pref=Pref.collapsed(["longitude"],iris.analysis.MEAN)

  P12 = P12.rolling_window("latitude",iris.analysis.MEAN,10)
  Pref = Pref.rolling_window("latitude",iris.analysis.MEAN,10)

  fig=plt.figure(figsize=(5,7))
  plt.subplot(311)
  plt.title("MC2")
  plt.xlim(-20,20)
  plt.plot([-15,-15,15,15,-15],[-1,16,16,-1,-1],"k",lw=1)
  plt.ylabel("Precip (mm/day)")
  plt.ylim(0,15)
  plt.subplot(312)
  plt.plot([-15,-15,15,15,-15],[-1,16,16,-1,-1],"k",lw=1)
  plt.title("MC12")
  plt.ylabel("Precip (mm/day)")
  plt.ylim(0,15)
  plt.xlim(-20,20)
  for i,year in enumerate(years):
    iplt.plot(P12.coord("latitude")[6:-6],P12[i,6:-6],lw=1)#,label="%04d-%02d"%(year,(year+1)%100))

  plt.subplot(313)
  plt.plot([-15,-15,15,15,-15],[-1,16,16,-1,-1],"k",lw=1)
  plt.xlim(-20,20)
  plt.ylim(0,15)
  plt.title("GPM-IMERG")
  plt.ylabel("Precip (mm/day)")
  for i,year in enumerate(years):
    iplt.plot(Pref.coord("latitude"),Pref[i],label="%04d-%02d"%(year,(year+1)%100),lw=1)

  fig.subplots_adjust(bottom=0.3,hspace=0.3)
  fig.legend(loc="lower center",ncol=4)
  fig.savefig(figname)
  plt.show()
