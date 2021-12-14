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
MC12_years = [2003,2014,2015,2016,2017]
ref_path = "/gws/nopw/j04/terramaris/emmah/era5/"
scratchpath = "/work/scratch-nopw/emmah/"
def load(year,MC):
    data=panMC(year,MC,"sea_surface").load_iris()
    data = data.extract("surface_net_downward_shortwave_flux")[0]+data.extract("surface_downward_heat_flux_in_sea_water_excludes_shortwave")[0]
    mean = data.collapsed("time",iris.analysis.MEAN)
    return mean

def load_ref(year,path):
    era5_12 = iris.load(path+"%04d_heat_rad_mon.nc"%(year),iris.Constraint(time=lambda t: t.point.month==12))
    era5_1  = iris.load(path+"%04d_heat_rad_mon.nc"%(year+1),iris.Constraint(time=lambda t: t.point.month==1))
    era5_2  = iris.load(path+"%04d_heat_rad_mon.nc"%(year+1),iris.Constraint(time=lambda t: t.point.month==2))
    for cube in era5_12+era5_1+era5_2:
      cube.coord("time").convert_units("days since 2003-01-01")
    iris.util.equalise_attributes(era5_12+era5_1+era5_2)
    era5 = (era5_12+era5_1+era5_2).merge()
    print(era5)
    idt_acc=iris.coords.AuxCoord(1/60/60/24,units="1/s")
    data = idt_acc*era5.extract("surface_net_downward_shortwave_flux")[0]\
          +idt_acc*era5.extract("surface_net_upward_longwave_flux")[0]\
          +idt_acc*era5.extract("surface_upward_latent_heat_flux")[0]\
          +idt_acc*era5.extract("surface_upward_sensible_heat_flux")[0]
    mean = data[:].collapsed("time",iris.analysis.MEAN)
    return mean 

def relax_callback(cube,field,filename):
  cube.var_name=None
  for coord in cube.coords():
     coord.var_name=None
 
def load_relax(year,path):
  relax=iris.load(["/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/%04d/KPPocean_relaxrun_inst_%04d%02d.nc"%x for x in [(year,year,12),(year,year+1,1),(year,year+1,2)]])
  relax = relax.concatenate()
  print(relax)
  data=relax.extract("surface_downward_heat_flux_in_sea_water_excludes_shortwave")[0]+relax.extract("surface_net_downward_shortwave_flux")[0]
  mean = data[:].collapsed("time",iris.analysis.MEAN)
  return mean
 
def main():
  MC12 = iris.cube.CubeList()
  ref = iris.cube.CubeList()
  relax = iris.cube.CubeList()

  for year in MC12_years:
    print(year)
    relax.append(load_relax(year,ref_path))
    ref.append(load_ref(year,ref_path))
    MC12.append(load(year,"MC12"))
    relax[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    ref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    MC12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")


  MC12=MC12.merge_cube().collapsed("time",iris.analysis.MEAN)

  ref=ref.merge_cube().collapsed("time",iris.analysis.MEAN)
  relax=relax.merge_cube().collapsed("time",iris.analysis.MEAN)

  MC12.data=MC12.data.filled(0)
  relax.data=relax.data.filled(0)

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

  ref = ref.regrid(MC12,iris.analysis.Linear())
  fig=bias_plots([MC12,relax,ref], ["MC12","Relax","Reference"] ,cmap1="bwr",cmap2="bwr",projection=ccrs.PlateCarree(),nmax=12,above=True,mark_inner=True)
  fig.set_figwidth(12)
  fig.set_figheight(7)
  fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
  fig.suptitle("Heat Flux into Ocean")
  plt.ion()
  plt.show()

if __name__=="__main__":
  main()
