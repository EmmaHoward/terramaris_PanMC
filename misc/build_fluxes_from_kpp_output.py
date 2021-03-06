import iris
import datetime as dt
import sys
from panMC import panMC
date = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))
datem1 =date+dt.timedelta(-6)

year = date.year - int(date.month < 6)

init = (date-dt.datetime(year,1,1)).days+1

dates = [date+dt.timedelta(i) for i in range(6)]
data = panMC(year,"MC12","sea_surface").load_iris(dates)
out = iris.cube.CubeList()

data.extract("surface_net_downward_shortwave_flux")[0].rename("swf")
data.extract("surface_downward_heat_flux_in_sea_water_excludes_shortwave")[0].rename("lwf")
data.extract("surface_downward_eastward_stress")[0].rename("taux")
data.extract("surface_downward_northward_stress")[0].rename("tauy")
data.extract("precipitation minus evaporation")[0].rename("pme")


out = data.extract(["swf","lwf","pme","taux","tauy"])

#out = iris.cube.CubeList([cube[23:-1] for cube in out])

for cube in out:
  cube.coords()[0].rename("time")
  cube.coord("time").points = cube.coord("time").points - 1/48

for name in ["shf","lhf"]:
  out.append(out.extract("lwf")[0].copy(data=0*out.extract("lwf")[0].data))
  out[-1].rename(name)

iris.save(out,"/work/scratch-nopw/emmah/%04d%02d%02d_surfacefluxes.nc"%(date.year,date.month,date.day))
