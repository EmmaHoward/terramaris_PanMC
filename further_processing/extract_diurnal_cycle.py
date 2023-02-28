#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python

#
# compute diurnal mean of precipitation
#

from panMC import panMC
from iris.experimental.equalise_cubes import equalise_attributes
import iris
from iris.coord_categorisation import add_categorised_coord
from subprocess import check_call
import sys

year=int(sys.argv[1])
MC = sys.argv[2]

path_template = "/gws/nopw/j04/terramaris/emmah/um/u-bs742/processed/MC_15_orog.nc"
scratchpath = "/work/scratch-pw2/emmah/tmp_%s/"%MC
finalpath = {"MC2":"/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/",
            "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/"}[MC]

check_call("mkdir -p %s"%scratchpath,shell=True)


def add_time_of_day(cube, coord, name='hour'):
    def _time_of_day(coord, value):
        pt = coord.units.num2date(value)
        return pt.hour+pt.minute/60.0+pt.second/3600
    add_categorised_coord(cube, name, coord,
                          _time_of_day)
						  
variables = ["convective_rainfall_amount","convective_snowfall_amount","stratiform_rainfall_amount","stratiform_snowfall_amount"]

def main(year,MC):
  data = panMC(year,MC,"rainfall").load_iris(variables=variables)
  for cube in data:
    add_time_of_day(cube,"time","hour")
  for month in [12,1,2]:
    ct = iris.Constraint(time=lambda t: t.point.month==month)
    new = iris.cube.CubeList([cube.extract(ct).aggregated_by("hour",iris.analysis.MEAN) for cube in data])
    iris.save(new,"%s/%04d%02d_diurnal_precip.nc"%(scratchpath,year+int(month<6),month),zlib=True)
    check_call("cp %s/%04d%02d_diurnal_precip.nc %s/postprocessed_outputs/precip/"%(scratchpath,year+int(month<6),month,finalpath),shell=True)
    check_call("rm %s/%04d%02d_diurnal_precip.nc"%(scratchpath,year+int(month<6),month),shell=True)

main(year,MC)


