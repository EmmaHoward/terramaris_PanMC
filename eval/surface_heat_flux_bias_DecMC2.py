import datetime as dt
import iris
from panMC import panMC
import cmocean
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import subprocess
from iris.coord_categorisation import add_hour,add_month
MC12_years = [2015,2016]
ref_path = "/gws/nopw/j04/terramaris/emmah/era5/"
scratchpath = "/work/scratch-nopw/emmah/"
relaxpath = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/"
MC2 = iris.cube.CubeList()
MC12 = iris.cube.CubeList()
ref = iris.cube.CubeList()
relax = iris.cube.CubeList()

def load(year,MC):
    dates = [dt.datetime(year,12,1)+dt.timedelta(i) for i in range(30)]
    data=panMC(year,MC,"sea_surface").load_iris(dates)
    data = data.extract("surface_net_downward_shortwave_flux")[0]+data.extract("surface_downward_heat_flux_in_sea_water_excludes_shortwave")[0]
    mean = data.collapsed("time",iris.analysis.MEAN)
    return mean

def load_ref(year,path):
    era5 = iris.load(path+"%04d_heat_rad_mon.nc"%(year),iris.Constraint(time=lambda t: t.point.month==12))
    print(era5)
    idt_acc=iris.coords.AuxCoord(1/60/60/24,units="1/s")
    data = idt_acc*era5.extract("surface_net_downward_shortwave_flux")[0]\
          +idt_acc*era5.extract("surface_net_upward_longwave_flux")[0]\
          +idt_acc*era5.extract("surface_upward_latent_heat_flux")[0]\
          +idt_acc*era5.extract("surface_upward_sensible_heat_flux")[0]
#    mean = data[:].collapsed("time",iris.analysis.MEAN)
    return data

def relax_callback(cube,field,filename):
  cube.var_name=None
  for coord in cube.coords():
     coord.var_name=None
 
def load_relax(year,path):
  relax=iris.load([path+"%04d/KPPocean_relaxrun_inst_%04d%02d.nc"%x for x in [(year,year,12)]])
  relax = relax.concatenate()
  print(relax)
  data=relax.extract("surface_downward_heat_flux_in_sea_water_excludes_shortwave")[0]+relax.extract("surface_net_downward_shortwave_flux")[0]
  mean = data[:].collapsed("time",iris.analysis.MEAN)
  return mean
 

for year in MC12_years:
  print(year)
  if year==2015:
    MC2.append(load(year,"MC2-tmp"))
  else:
    MC2.append(load(year,"MC2"))
  relax.append(load_relax(year,relaxpath))
  ref.append(load_ref(year,ref_path))
  MC12.append(load(year,"MC12"))
  relax[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  ref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  MC12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  MC2[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")


MC2=MC2.merge_cube().collapsed("time",iris.analysis.MEAN)
MC12=MC12.merge_cube().collapsed("time",iris.analysis.MEAN)
iris.util.equalise_attributes(ref)
ref=ref.merge_cube().collapsed("time",iris.analysis.MEAN)
relax=relax.merge_cube().collapsed("time",iris.analysis.MEAN)

MC12.data=MC12.data.filled(0)
relax.data=relax.data.filled(0)
MC2.data=MC2.data.filled(0)

ref = ref.regrid(MC12,iris.analysis.Linear())
MC2 = MC2.regrid(MC12,iris.analysis.Linear(extrapolation_mode="mask"))

fig=bias_plots([MC2,MC12,relax,ref], ["MC2","MC12","Relax","Reference"] ,cmap1="bwr",cmap2="bwr",projection=ccrs.PlateCarree(),nmax=12,above=True,mark_inner=True,minmaxN1rangeN2=(-130,130,14,130,14))
fig.set_figwidth(12)
fig.set_figheight(7)
fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
fig.suptitle("Heat Flux into Ocean")
plt.ion()
plt.show()

