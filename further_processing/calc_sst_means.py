#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
import iris
import iris.cube
from panMC import panMC
import sys
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import iris.plot as iplt
import subprocess
from iris.coord_categorisation import add_hour,add_month,add_day_of_year
import datetime as dt

job = int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
year = int(sys.argv[1])
MC = sys.argv[2]

MC12_years = [2003,2014,2015,2016,2017,2018]
ref_path = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/"
scratchpath = "/work/scratch-pw2/emmah/tmp_%s/"%MC
sst12 = iris.cube.CubeList()
sstref = iris.cube.CubeList()
finalpath = {"MC2": "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/sst/",
             "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/sst/"}[MC]

cx = iris.Constraint(longitude=lambda lon: 90<=lon<=155)
cy = iris.Constraint(latitude=lambda lat: -15<=lat<=15)

def adjust_doyr(cube):
  doyr = cube.coord("doyr").points
  if (cube.coord("year").points%4==0).all():
    doyr[doyr>180] -= 1
  doyr[doyr<180] += 365
  cube.coord("doyr").points = doyr

def calc_diurnal_mean_sst(year,MC,path):
    outname = "%s_%04d%02d_diurnal_sea_temperature.nc"%(MC,year,(year+1)%100)
    data=panMC(year,MC,"sea_water_temperature").load_iris()[0]
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    data.coord("time").bounds = np.array([np.round(data.coord("time").points)-1,np.round(data.coord("time").points)]).T
    data.coord("time").points = np.round(data.coord("time").points)-0.5
    add_hour(data,"time","hour")
    add_month(data,"time","month")
    dcmean = data[:,:6].aggregated_by(["hour","month"],iris.analysis.MEAN)
    iris.save(dcmean,scratchpath+outname,zlib=True)
    subprocess.check_call("cp %s%s %s%s"%(scratchpath,outname,path,outname),shell=True)
    subprocess.check_call("rm %s%s"%(scratchpath,outname),shell=True)

def calc_daily_mean_sst(year,MC,path):
    outname = "%s_%04d%02d_daily_sea_temperature.nc"%(MC,year,(year+1)%100)
    data=panMC(year,MC,"sea_water_temperature").load_iris()[0]
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    data.coord("time").bounds = np.array([np.round(data.coord("time").points)-1,np.round(data.coord("time").points)]).T
    data.coord("time").points = np.round(data.coord("time").points)-0.5
    add_day_of_year(data,"time","doyr")
    dmean = data[:,:6].aggregated_by(["doyr"],iris.analysis.MEAN)
    iris.save(dmean,scratchpath+outname,zlib=True)
    subprocess.check_call("cp %s%s %s%s"%(scratchpath,outname,path,outname),shell=True)
    subprocess.check_call("rm %s%s"%(scratchpath,outname),shell=True)
 
def calc_diurnal_range_sst(year,MC,path):
    outname = "%s_%04d%02d_sst_diurnal_range.nc"%(MC,year,(year+1)%100)
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


if job == 0:
  calc_diurnal_mean_sst(year,MC,finalpath)
elif job == 1:
  calc_diurnal_range_sst(year,MC,finalpath)
elif job == 2:
  calc_daily_mean_sst(year,MC,finalpath)
