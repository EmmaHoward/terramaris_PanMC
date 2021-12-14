import iris
import datetime as dt
import sys
date = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))
datem1 =date+dt.timedelta(-6)

year = date.year - int(date.month < 6)

init = (date-dt.datetime(year,1,1)).days+1

path = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/u-cf309/%04d/%04d%02d%02dT0000Z/"%(year,date.year,date.month,date.day)
#path2 = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/u-cf309/%04d/%04d%02d%02dT0000Z/"%(year,datem1.year,datem1.month,datem1.day)
data = iris.load([path+"KPPocean_?????.nc"])#,path2+"KPPocean_%05d.nc"%(init-1)])
data = data.concatenate()

out = iris.cube.CubeList()

data.extract("solar_in")[0].rename("swf")
data.extract("nsolar_in")[0].rename("lwf")
data.extract("taux_in")[0].rename("taux")
data.extract("tauy_in")[0].rename("tauy")
data.extract("PminusE_in")[0].rename("pme")


out = data.extract(["swf","lwf","pme","taux","tauy"])

#out = iris.cube.CubeList([cube[23:-1] for cube in out])

for cube in out:
  cube.coords()[0].rename("time")
  cube.coord("time").convert_units("days since %04d-01-01 00:00:00"%year)
  cube.coord("time").points = cube.coord("time").points - 1/48

for name in ["shf","lhf"]:
  out.append(out.extract("lwf")[0].copy(data=0*out.extract("lwf")[0].data))
  out[-1].rename(name)

iris.save(out,"/work/scratch-nopw/emmah/%04d%02d%02d_surfacefluxes.nc"%(date.year,date.month,date.day))
