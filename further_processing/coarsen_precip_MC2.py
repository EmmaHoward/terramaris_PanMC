#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python

from panMC import panMC
import iris
from iris.coord_categorisation import add_hour,add_day_of_year
import datetime as dt
from subprocess import check_call
import sys
import os

scratchpath = "/work/scratch-pw2/emmah/tmp_MC2/"
path_template = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"
finalpath =  "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"

def coarsen(year,month):
  t0 = dt.datetime(year,12,1)
  dates = [t0 + dt.timedelta(i) for i in range(90)]
  dates = [t for t in dates if t.month==month]
  data = panMC(year,"MC2","rainfall").load_iris(dates,variables=["stratiform_rainfall_amount"])
  print(data)
  data=data[0]
  add_hour(data,"time","hour")
  add_day_of_year(data,"time","doyr")
  agg = data.aggregated_by(["hour","doyr"],iris.analysis.MEAN)
  template = iris.load(path_template)[0][54:-55,42:-42] # bounds extract inner domain
  template.coord("longitude").guess_bounds()
  template.coord("latitude").guess_bounds()
  agg.coord("longitude").guess_bounds()
  agg.coord("latitude").guess_bounds()
  regrid = agg.regrid(template,iris.analysis.AreaWeighted())
  filename = "%04d%02d_hourly_coarsened_rainfall.nc"%(year+int(month<6),month)
  iris.save(regrid,"%s/%s"%(scratchpath,filename))
  check_call("cp %s/%s %s"%(scratchpath,filename,finalpath),shell=True)
  check_call("rm %s/%s"%(scratchpath,filename),shell=True)
 
year=int(sys.argv[1])
job =int(os.environ["SLURM_ARRAY_TASK_ID"]) -1 
month = [12,1,2][job]
coarsen(year,month)
