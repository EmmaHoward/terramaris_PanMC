#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python

#
# regrid daily mean 850 and 200 hPa winds to era5 grid
#

from matplotlib.colors import ListedColormap
from panMC import panMC
import seaborn as sns
import iris
from scipy.fftpack import dct
import cftime
from iris.experimental.equalise_cubes import equalise_attributes
#from panMC import panMC
from iris.coord_categorisation import add_day_of_year
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import iris.cube
import sys
import os
from subprocess import check_call
cy =iris.Constraint(latitude=lambda y: -15<=y<=15)
cx =iris.Constraint(longitude=lambda y: 90<=y<=155)
cz = iris.Constraint(pressure=lambda p: p in [200,850])

year = int(sys.argv[1])
MC = sys.argv[2]
scratchpath = "/work/scratch-pw2/emmah/tmp_%s/"%MC
finalpath = {"MC2": "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/wind/",
             "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/wind/"}[MC]


def adjust_doyr(cube):
  doyr = cube.coord("doyr").points
  if (cube.coord("year").points%4==0).all():
    doyr[doyr>180] -= 1
  doyr[doyr<180] += 365
  cube.coord("doyr").points = doyr

    
def regrid_wind(year,MC):
    print(year)
    yc = iris.coords.DimCoord(year,long_name="year",units="years")
    u = panMC(year,MC,"eastward_wind_pd").load_iris(Constraints=cz)[0]
    v = panMC(year,MC,"northward_wind_pd").load_iris(Constraints=cz)[0]
    u.coord("time").points = u.coord("time").points - 1/48
    v.coord("time").points = v.coord("time").points - 1/48
    add_day_of_year(u,"time","doyr")
    add_day_of_year(v,"time","doyr")
    u.coord("time").points = u.coord("time").points + 1/48
    v.coord("time").points = v.coord("time").points + 1/48
    u = u.aggregated_by("doyr",iris.analysis.MEAN)
    v = v.aggregated_by("doyr",iris.analysis.MEAN)
#
    template =  iris.load("/gws/nopw/j04/terramaris/emmah/era5/uv_201415.nc",cx&cy).extract("eastward_wind")[0]
    template = template.rolling_window("latitude",iris.analysis.MEAN,2)
    template.coord("longitude").guess_bounds()
#
    u.coord("longitude").guess_bounds()
    u.coord("latitude").guess_bounds()
    v.coord("longitude").guess_bounds()
    v.coord("latitude").guess_bounds()
    u.coord("longitude").coord_system = template.coord("longitude").coord_system
    u.coord("latitude").coord_system  = template.coord("latitude").coord_system
    v.coord("longitude").coord_system = template.coord("longitude").coord_system
    v.coord("latitude").coord_system  = template.coord("latitude").coord_system
#
    u = u.regrid(template,iris.analysis.AreaWeighted())
    v = v.regrid(template,iris.analysis.AreaWeighted())
    u.add_aux_coord(yc)
    v.add_aux_coord(yc)
    adjust_doyr(u)
    adjust_doyr(v)
#
    iris.save([u,v],"%s/%s_%04d%02d_winds_regridded0p25.nc"%(scratchpath,MC,year,(year+1)%100),zlib=True)
    check_call("cp %s/%s_%04d%02d_winds_regridded0p25.nc %s"%(scratchpath,MC,year,(year+1)%100,finalpath),shell=True) 
    check_call("rm %s/%s_%04d%02d_winds_regridded0p25.nc"%(scratchpath,MC,year,(year+1)%100),shell=True)

regrid_wind(year,MC)
