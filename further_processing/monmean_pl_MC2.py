#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python

# 
# compute monthly means of MC2 pressure level data, coarsened to MC12 grid
#

import matplotlib.pyplot as plt
import os
import numpy as np
import iris
import cmocean
from iris.coord_categorisation import add_month
import sys
from panMC import panMC
from subprocess import check_call
import iris.cube
import datetime as dt
year=int(sys.argv[1])
MC = sys.argv[2]
task = sys.argv[3]

finalpath = "/gws/nopw/j04/terramaris/panMC_um/%s_*/postprocessed_outputs/monmean_pl/"%MC
path_template = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"
scratchpath = "/work/scratch-pw2/emmah/tmp_%s/"%MC

def average(scratchpath, year, mon, MC, var,j):
    dates = [dt.datetime(year+int(mon<6),mon,1)+ dt.timedelta(i) for i in range({12:31,1:31,2:28}[mon])]
    H = panMC(year,MC,"Heavyside").load_iris(dates)[0][:,j]
    template = iris.load(path_template)[0][54:-55,42:-42] # bounds extract inner domain
    template.coord("latitude").guess_bounds()
    template.coord("longitude").guess_bounds()
    data=panMC(year,MC,var).load_iris(dates)[0][:,j]
    if MC == "MC2":
      data.data.mask += (H.data<1)
      data = data.collapsed("time",iris.analysis.MEAN)
      data.coord("longitude").guess_bounds()
      data.coord("latitude").guess_bounds()
      data = data.regrid(template,iris.analysis.AreaWeighted())
    else:
      data = data.collapsed("time",iris.analysis.MEAN)
    iris.save(data,"%s/%s_pressurelevel_%d%02d_%s_%d.nc"%(scratchpath,MC,year,mon,var,j),zlib=True)
   

def main(scratchpath, mon, year, MC, variables,j):
  for var in variables:
    for mon in [12,1,2]:
      average(scratchpath, year, mon, MC, var,j)


def collate(scratchpath, year, MC, variables):
  files  = []
  for mon in [12,1,2]:
    for var in variables:
      for j in range(2,22):
         files.append("%s/%s_pressurelevel_%d%02d_%s_%d.nc"%(scratchpath,MC,year,mon,var,j))
  data = iris.load(files)
  a = iris.cube.CubeList()
  for cube in data:
    if cube.ndim==4:
      for i in range(cube.shape[0]):
        a.append(cube[i])
    else:
      a.append(cube)
#  import pdb;pdb.set_trace()
  data=a.merge()
  iris.save(data,"%s/%s_monmean_pl_%d.nc"%(scratchpath,MC,year),zlib=True)
  check_call("cp %s/%s_monmean_pl_%d.nc %s"%(scratchpath,MC,year,finalpath),shell=True)
  check_call("rm %s/%s_pressurelevel_%d*.nc"%(scratchpath,MC,year),shell=True)


if task=="calc":
  job = int(os.environ["SLURM_ARRAY_TASK_ID"]) -1
  j = job%22
  mon=job//22
  main(scratchpath,mon, year,MC,['omega', 'specific_humidity_pd',  'theta_w','northward_wind_pd','temperature_pd','eastward_wind_pd'],j)
elif task=="collate":
  collate(scratchpath,year,MC,['omega', 'specific_humidity_pd', 'eastward_wind_pd', 'theta_w','northward_wind_pd','temperature_pd'])
else:
  print("3rd arg should be calc or collate")
  exit()
