#
# compute monthly means of MC12 pressure level data
#

#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
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

finalpath = "/gws/nopw/j04/terramaris/panMC_um/%s_*/postprocessed_outputs/monmean_pl/"%MC
path_template = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"
scratchpath = "/work/scratch-pw2/emmah/tmp_%s/"%MC

def average(scratchpath, year, mon, MC, var):
    dates = [dt.datetime(year+int(mon<6),mon,1)+ dt.timedelta(i) for i in range({12:31,1:31,2:28}[mon])]
    if MC=="MC2":
        H = panMC(year,MC,"Heavyside").load_iris(dates)[0]
        template = iris.load(path_template)[0][54:-55,42:-42] # bounds extract inner domain
        template.coord("latitude").guess_bounds()
        template.coord("longitude").guess_bounds()

    data=panMC(year,MC,var).load_iris(dates)[0]
    if MC == "MC2":
      if var in ['eastward_wind_pd', 'northward_wind_pd']:
        H2 = H.regrid(data,iris.analysis.Linear)
      else:
        H2 = H.copy()
      data.data.mask += (H2.data>0)
      data = data.collapsed("time",iris.analysis.MEAN)
      data.coord("longitude").guess_bounds()
      data.coord("latitude").guess_bounds()
      data = data.regrid(template,iris.analysis.AreaWeighted())
    else:
      data = data.collapsed("time",iris.analysis.MEAN)
    iris.save(data,"%s/%s_pressurelevel_%d%02d_%s.nc"%(scratchpath,MC,year,mon,var),zlib=True)
   

def main(scratchpath, year, MC, variables):
  for var in variables:
    for mon in [12,1,2][:]:
      average(scratchpath, year, mon, MC, var)
  data = iris.load(["%s/%s_pressurelevel_%d*_%s.nc"%(scratchpath,MC,year,var) for var in variables])
  data=data.merge()
  iris.save(data,"%s/%s_monmean_pl_%d.nc"%(scratchpath,MC,year),zlib=True)
  check_call("cp %s/%s_monmean_pl_%d.nc %s"%(scratchpath,MC,year,finalpath),shell=True)
#  check_call("rm %s/%s_monmean_pl_%d*.nc"%(scratchpath,MC,year),shell=True)

main(scratchpath,year,MC,['temperature_pd','omega', 'specific_humidity_pd', 'theta_w', 'eastward_wind_pd', 'northward_wind_pd'][:])

