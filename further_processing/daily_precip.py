#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python

#
# Compute daily mean rainfall 
#

from panMC import panMC
from iris.experimental.equalise_cubes import equalise_attributes
import iris
from iris.coord_categorisation import add_day_of_year
from subprocess import check_call
import sys

year=int(sys.argv[1])
MC = sys.argv[2]

path_template = "/gws/nopw/j04/terramaris/emmah/um/u-bs742/processed/MC_15_orog.nc"
scratchpath = "/work/scratch-pw2/emmah/tmp_%s/"%MC
finalpath = {"MC2":"/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/",
            "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/"}[MC]

check_call("mkdir -p %s"%scratchpath,shell=True)


variables = ["convective_rainfall_amount","convective_snowfall_amount","stratiform_rainfall_amount","stratiform_snowfall_amount"]

def main(year,MC):
  data = panMC(year,MC,"rainfall").load_iris(variables=variables)
  if MC == "MC2":
    cube = data[0]+data[1]
  else:
    cube = data[0]+data[1]+data[2]+data[3]
  add_day_of_year(cube,"time","doyr")
  cube = cube.aggregated_by("doyr",iris.analysis.SUM)  
  iris.save(cube,"%s/%04d%02d_daily_precip.nc"%(scratchpath,year,(year+1)%100),zlib=True)
  check_call("cp %s/%04d%02d_daily_precip.nc %s/postprocessed_outputs/precip/"%(scratchpath,year,(year+1)%100,finalpath),shell=True)
  check_call("rm %s/%04d%02d_daily_precip.nc"%(scratchpath,year,(year+1)%100),shell=True)

main(year,MC)


