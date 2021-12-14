import iris
from panMC import panMC
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import os
import numpy as np
import datetime as dt
import subprocess
from iris.coord_categorisation import add_hour,add_month
MC12_years = [2003,2014,2015,2016,2017,2018]
MC2_years = [2016]
ref_path = "/gws/nopw/j04/terramaris/emmah/monthly_cmems/hmix/"
scratchpath = "/work/scratch-nopw/emmah/"
MC12 = iris.cube.CubeList()
MC2 = iris.cube.CubeList()
ref = iris.cube.CubeList()


def load(year,MC):
   data=panMC(year,MC,"mixedlayerdepth").load_iris()[0]
   print(data)
   dcmean = data.collapsed("time",iris.analysis.MEAN)
   return dcmean

def load_ref(year,path):
    data = iris.load_cube(path+"GLOBAL_REANALYSIS_PHY_001_030-TDS_%04d-12-16.nc"%(year))
    mean = data.collapsed("time",iris.analysis.MEAN)
    mean.coord("time").attributes = {}  
    return mean 

for year in MC2_years:
  print(year)
  ref.append(load_ref(year,ref_path))
  MC12.append(load(year,"MC12"))
  if year==2015:
    MC2.append(load(year,"MC2-tmp"))
  else:
    MC2.append(load(year,"MC2"))
  ref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  MC12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  MC2[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")


MC12=MC12.merge_cube()#.collapsed("time",iris.analysis.MEAN)
MC2=MC2.merge_cube()#.collapsed("time",iris.analysis.MEAN)

ref=ref.merge_cube()#.collapsed("time",iris.analysis.MEAN)

ref.coord("longitude").var_name=None
MC12.coord("longitude").long_name=None
MC12.coord("longitude").attributes={}
ref.coord("longitude").attributes={}
MC12.coord("longitude").guess_bounds()

ref.coord("latitude").var_name=None
MC12.coord("latitude").long_name=None
MC12.coord("latitude").attributes={}
ref.coord("latitude").attributes={}
MC12.coord("latitude").guess_bounds()

ref = ref.regrid(MC2,iris.analysis.Linear())
MC12 = MC12.regrid(MC2,iris.analysis.Nearest())
#relax = relax.regrid(MC2,iris.analysis.Nearest())

fig=bias_plots([MC2,MC12,ref], ["MC2","MC12","Reference"] ,cmap1="viridis",cmap2="bwr",projection=ccrs.PlateCarree(),nmax=12,above=True,mark_inner=True)
fig.set_figwidth(12)
fig.set_figheight(7)
fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
fig.suptitle("Surface Heat Flux Imbalance")
plt.ion()
plt.show()

