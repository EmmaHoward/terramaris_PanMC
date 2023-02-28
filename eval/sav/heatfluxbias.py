#
# supplementary figure: biases in heat fluxes into ocean
# compares MC2, MC12, relaxation runs and observations
#

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
ref_path = "/gws/nopw/j04/terramaris/emmah/era5/"
relaxpath = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/"
relaxyears =range(2003,2018)
def load(year,MC):
    # load heat fluxes from KPP output (which splits surface forcing into shortwave and everything else)
    data=panMC(year,MC,"sea_surface").load_iris()
    data = data.extract("surface_net_downward_shortwave_flux")[0]+data.extract("surface_downward_heat_flux_in_sea_water_excludes_shortwave")[0]
    mean = data.collapsed("time",iris.analysis.MEAN)
    return mean

def load_ref(year,path):
    # load era5 heat fluxes a
    era5 = iris.load(path+"%04d_heat_rad_mon.nc"%(year))
    print(era5)
    idt_acc=iris.coords.AuxCoord(1/60/60/24,units="1/s") # converts units to W/m2
    data = idt_acc*era5.extract("surface_net_downward_shortwave_flux")[0]\
          +idt_acc*era5.extract("surface_net_upward_longwave_flux")[0]\
          +idt_acc*era5.extract("surface_upward_latent_heat_flux")[0]\
          +idt_acc*era5.extract("surface_upward_sensible_heat_flux")[0]
    mean = data[:].collapsed("time",iris.analysis.MEAN)
    return mean

def relax_callback(cube,field,filename):
  # loading metadata helper function
  cube.var_name=None
  for coord in cube.coords():
     coord.var_name=None
 
def load_relax(year,path):
  # load fluxes from relaxation runs
  relax=iris.load([path+"%04d/KPPocean_relaxrun_inst_%04d%02d.nc"%x for x in [(year,year,12),(year,year+1,1),(year,year+1,2)]])
  relax = relax.concatenate()
  print(relax)
  data=relax.extract("surface_downward_heat_flux_in_sea_water_excludes_shortwave")[0]+relax.extract("surface_net_downward_shortwave_flux")[0]
  mean = data[:].collapsed("time",iris.analysis.MEAN)
  return mean
 
def main(MC12_years,MC2_years,figname):
  #load data and make plots
  MC2 = iris.cube.CubeList()
  MC12 = iris.cube.CubeList()
  ref = iris.cube.CubeList()
  relax = iris.cube.CubeList()
  for year in MC12_years:
    print(year)
    ref.append(load_ref(year,ref_path))
    ref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    MC12.append(load(year,"MC12"))
    MC12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  for year in relaxyears:
    relax.append(load_relax(year,relaxpath))
    relax[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  for year in MC2_years:
    MC2.append(load(year,"MC2"))
    MC2[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  if len(MC12_years)==1:
    MC12=MC12.merge_cube()#.collapsed("time",iris.analysis.MEAN)
    iris.util.equalise_attributes(ref)
    ref=ref.merge_cube()#.collapsed("time",iris.analysis.MEAN)
    relax=relax.merge_cube()#.collapsed("time",iris.analysis.MEAN)
  else:
    MC12=MC12.merge_cube().collapsed("time",iris.analysis.MEAN)
    iris.util.equalise_attributes(ref)
    ref=ref.merge_cube().collapsed("time",iris.analysis.MEAN)
    relax=relax.merge_cube().collapsed("time",iris.analysis.MEAN)
  if len(MC2_years)==1:
    MC2=MC2.merge_cube()#.collapsed("time",iris.analysis.MEAN)
  else:
    MC2=MC2.merge_cube().collapsed("time",iris.analysis.MEAN)
  MC12.data=MC12.data.filled(0)
  relax.data=relax.data.filled(0)
  MC2.data=MC2.data.filled(0)

  ref = ref.regrid(MC12,iris.analysis.Linear())
  MC2 = MC2.regrid(MC12,iris.analysis.Linear(extrapolation_mode="mask"))

  fig=bias_plots([MC2,MC12,relax,ref], ["MC2","MC12","Relax","Reference"] ,cmap1="bwr",cmap2="BrBG_r",projection=ccrs.PlateCarree(),nmax=12,above=True,mark_inner=True,minmaxN1rangeN2=(-130,130,14,65,14))
  fig.set_figwidth(12)
  fig.set_figheight(7)
  fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
  fig.suptitle("Heat Flux into Ocean")
  fig.savefig(figname)

