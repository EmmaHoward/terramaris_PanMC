import iris
from panMC import panMC
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import subprocess
from iris.coord_categorisation import add_day_of_year
MC12_years = [2015]#[2003,2014,2015,2016,2017,2018]
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monthly_diurnal_sst/"
ref_path = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/"
scratchpath = "/work/scratch-pw2/emmah/"
sst12 = iris.cube.CubeList()
sstref = iris.cube.CubeList()

def load(year,MC,path):
  outname = "%s_%04d%02d_sst_diurnal_range.nc"%(MC,year,(year+1)%100)
  print(year)
  if 0:#os.path.exists(path+outname):
    print("exists")
    dcmean = iris.load(path+outname)[0]
  else:
    data=panMC(year,MC,"sea_water_temperature").load_iris()[0]
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    data.coord("time").bounds = np.array([np.round(data.coord("time").points)-1,np.round(data.coord("time").points)]).T
    data.coord("time").points = np.round(data.coord("time").points)-0.5
    add_day_of_year(data,"time","doyr")
    dcmax = data[:,0].aggregated_by("doyr",iris.analysis.MAX)
    dcmean = data[:,3].aggregated_by("doyr",iris.analysis.MEAN)
    dcmin = data[:,0].aggregated_by("doyr",iris.analysis.MIN)
    dc = dcmax - dcmin
    dc.rename("sst_diurnal_range")
    iris.save([dc,dcmean],scratchpath+outname,zlib=True)
    subprocess.check_call("cp %s%s %s%s"%(scratchpath,outname,path,outname),shell=True)
    subprocess.check_call("rm %s%s"%(scratchpath,outname),shell=True)
  return dcmean

for year in MC12_years:
  load(year,"MC12",MC12_path)
