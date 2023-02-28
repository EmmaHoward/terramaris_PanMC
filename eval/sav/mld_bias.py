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
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monthly_diurnal_sst/"
ref_path = "/gws/nopw/j04/terramaris/emmah/monthly_cmems/hmix/"
scratchpath = "/work/scratch-nopw/emmah/"
def load(year,MC):
    data=panMC(year,MC,"mixedlayerdepth").load_iris()[0]
    dcmean = data.collapsed("time",iris.analysis.MEAN)
    return dcmean

def load_ref(year,path):
    data = iris.load_cube(path+"GLOBAL_REANALYSIS_PHY_001_030-TDS_%04d-12-16.nc"%(year))
    mean = data[:].collapsed("time",iris.analysis.MEAN)
    mean.coord("time").attributes = {}  
    return mean 

def main(MC12_years,MC2_years,scratchpath,figname=None):
  MC2 = iris.cube.CubeList()
  MC12 = iris.cube.CubeList()
  ref = iris.cube.CubeList()
  for year in MC12_years:
    print(year)
    ref.append(load_ref(year,ref_path))
    MC12.append(load(year,"MC12"))
    ref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    MC12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  for year in MC2_years:
    MC2.append(load(year,"MC2"))
    MC2[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  MC12=MC12.merge_cube().collapsed("time",iris.analysis.MEAN)
  ref=ref.merge_cube().collapsed("time",iris.analysis.MEAN)
  MC2=MC2.merge_cube().collapsed("time",iris.analysis.MEAN)

  ref.coord("longitude").var_name=None
  MC12.coord("longitude").long_name=None
  MC2.coord("longitude").long_name=None
  MC12.coord("longitude").attributes={}
  MC2.coord("longitude").attributes={}
  ref.coord("longitude").attributes={}
  MC12.coord("longitude").guess_bounds()
  MC2.coord("longitude").guess_bounds()

  ref.coord("latitude").var_name=None
  MC2.coord("latitude").long_name=None
  MC12.coord("latitude").long_name=None
  MC2.coord("latitude").attributes={}
  MC12.coord("latitude").attributes={}
  ref.coord("latitude").attributes={}
  MC2.coord("latitude").guess_bounds()
  MC12.coord("latitude").guess_bounds()
  ref.rename("cmems hmix")
  MC12.rename("MC12 hmix")
  MC2.rename("MC2 hmix")
  iris.save([MC12,ref],"%s/mld.nc"%scratchpath,zlib=True)
  ref = ref.regrid(MC12,iris.analysis.Linear())
  MC2 = MC2.regrid(MC12,iris.analysis.Linear(extrapolation_mode="mask"))
  fig=bias_plots([MC2,MC12,ref], ["MC2","MC12","Reference"] ,cmap1="cmo.deep",cmap2="PiYG",projection=ccrs.PlateCarree(),nmax=12,above=True,mark_inner=True,minmaxN1rangeN2=(0,70,8,22.5,10))
  fig.set_figwidth(12)
  fig.set_figheight(7)
  fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
  fig.tight_layout()
  fig.suptitle("Mixed Layer Depth")
  fig.savefig(figname)

