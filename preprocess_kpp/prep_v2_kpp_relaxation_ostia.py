#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#
# Generate 30 day rolling mean of ocean reanalysis on KPP grid 
#
import sys
from gridfill import fill_cube
import os
from progressbar import progressbar
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np
import scipy.ndimage as nd
import datetime as dt
# Any data on the KPP grid
def prep(job,year1):
  year2=year1+1
  print(job)
  template = iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_initcond/20161101_initcond_N1280_kmax.nc")
  t0 = dt.datetime(year1,10,27)+dt.timedelta(job)
  t_read = dt.datetime(year1,10,27)+dt.timedelta(5*(job//5))
  ct = iris.Constraint(time = lambda t: t.point.day==t0.day and t.point.month==t0.month)
#job=int(sys.argv[1])
# Input data
#data = iris.load(["/work/n02/n02/emmah/cmems/201516/GLOBAL_REANALYSIS_PHY_001_030-TDS_2015-10-*.nc","/work/n02/n02/emmah/cmems/201516/GLOBAL_REANALYSIS_PHY_001_030-TDS_2015-11-*.nc"],["sea_water_salinity","sea_water_potential_temperature"])
  data2  = iris.load("/work/scratch-nopw/emmah/cmems/%04d%02d/GLOBAL_REANALYSIS_PHY_001_030-TDS_%04d-%02d-%02d.nc"%(year1,year2%100,t_read.year,t_read.month,t_read.day),ct).extract(["sea_water_salinity","sea_water_potential_temperature"])
# Artifact of how I downloaded cmems
#data2=data
  equalise_attributes(data2)
# remove attributes specific to cmems file structure, flip depth coord
  for j,cube in enumerate(data2):
    del cube.coord('time').attributes["valid_min"]
    del cube.coord('time').attributes["valid_max"]
    del cube.coord('depth').attributes["valid_min"]
    del cube.coord('depth').attributes["valid_max"]
    del cube.coord('depth').attributes["positive"]
    del cube.coord('depth').attributes["_CoordinateZisPositive"]
    cube.coord("depth").points = -1*cube.coord("depth").points
    cube.data=np.ma.masked_array(cube.data,mask=cube.data.mask+np.isnan(cube.data))
  for cube in data2:
    cube.coord("longitude").guess_bounds()
    cube.coord("latitude").guess_bounds()
  print(data2)
# merge
# Take 30 day running mean - note this discards first and last 15 days, so download extra
  data2 = iris.cube.CubeList([cube.regrid(template[0],iris.analysis.AreaWeighted()) for cube in data2])

# fill vertically with lowest valid point
# then, replace top 5 m with 5m value


# Vertical coordinate (NZP1)
  z = np.array([-0.6390094, -1.925694, -3.229947, -4.552257, -5.893128, -7.253092,
    -8.632697, -10.03252, -11.45317, -12.89526, -14.35946, -15.84645,
    -17.35696, -18.89175, -20.45159, -22.03733, -23.64985, -25.29005,
    -26.95892, -28.65745, -30.38675, -32.14793, -33.94219, -35.7708,
    -37.63511, -39.53653, -41.47657, -43.45684, -45.47903, -47.54497,
    -49.65658, -51.81593, -54.02523, -56.28684, -58.60331, -60.97736,
    -63.41193, -65.91018, -68.47553, -71.1117, -73.82271, -76.61291,
    -79.48708, -82.45045, -85.50872, -88.66817, -91.93572, -95.31905,
    -98.82664, -102.468, -106.2537, -110.1958, -114.3075, -118.6044,
    -123.1037, -127.8256, -132.7932, -138.0335, -143.5783, -149.465,
    -155.7388, -162.4542, -169.6781, -177.4939, -186.0077, -195.3565,
    -205.7225, -217.3551, -230.6084, -244.7767])

# Interpolate to vertical coordinate
  sample_points = [('depth', z)]
  data2 = iris.cube.CubeList([cube.interpolate(sample_points, iris.analysis.Linear()) for cube in data2])

  z_3d = (np.ones((70,428,538)).T*(z)).T

  bathy = iris.load("lsm_ocndepth_kmax.nc","max_depth")[0]

  data2[0].data = data2[0].data.filled(0)*(z_3d >= bathy.data)
  data2[1].data = data2[1].data.filled(0)*(z_3d >= bathy.data)
  indices = iris.load("kmax_indices.nc")
  ix = indices.extract("x-indices")[0].data
  iy = indices.extract("y-indices")[0].data

  d1 =data2[0].data
  d2 =data2[1].data

  d1[0] = d1[3]
  d1[1] = d1[3]
  d1[2] = d1[3]

  d2[0] =d2[3]
  d2[1] =d2[3]
  d2[2] =d2[3]

  ny,nx=bathy.shape
  for i in progressbar(range(ny)):
    for j in range(nx):
      if ix[i,j]>0:  
        d1[:,i,j] = d1[:,int(iy[i,j]),int(ix[i,j])]
        d2[:,i,j] = d2[:,int(iy[i,j]),int(ix[i,j])]
 
  data2[0].data=d1
  data2[1].data=d2

#x=np.ma.masked_array(data2[0][0].data,mask=(z_3d <= bathy.data))
#print(((x==0).any(axis=0)).sum())
#print(((iy>0)).sum())

#import matplotlib.pyplot as plt
#plt.pcolormesh((x==0).sum(axis=0)*(iy>0))
#plt.show()

#exit()

# Change names and units to match KPP requirements
  s = data2.extract("sea_water_salinity")[0]
  s.rename("salinity")
  s.coord('depth').rename("z")
  s.coord('time').convert_units("days since %04d-01-01 00:00:00"%year1)
  s.coord('time').rename('t')
  #s.coord('t').points = s.coord('t').points+0.5
  t = data2.extract("sea_water_potential_temperature")[0]
  t.rename("temperature")
  t.coord('depth').rename("z")
  t.coord('time').convert_units("days since %04d-01-01 00:00:00"%year1)
  t.coord('time').rename('t')
#t.coord('t').points = t.coord('t').points+0.5
  iris.save([t,s],"/work/scratch-nopw/emmah/cmems/%04d%02d/%04d%02d%02d_relax_kmax.nc"%(year1,year2%100,t0.year,t0.month,t0.day))
#iris.save(t,"/gws/nopw/j04/terramaris/emmah/cmems_ostia/%04d11_temperature.nc"%year1)



def roll(job,year1):
  year2 = year1+1
  t0 = dt.datetime(year1,11,1)+dt.timedelta(job)
  if t0.month==3 and t0.day > 17:
    return
  times =  [t0+dt.timedelta(i) for i in range(-3,4)]
  os.chdir("/work/scratch-nopw/emmah/cmems/")
  data = iris.load(["/work/scratch-nopw/emmah/cmems/%04d%02d/%04d%02d%02d_relax_kmax.nc"%(year1,year2%100,t1.year,t1.month,t1.day) for t1 in times])
  data = data.merge()
  t= data.extract("temperature")[0].collapsed("t",iris.analysis.MEAN)
  s= data.extract("salinity")[0].collapsed("t",iris.analysis.MEAN)
#  iris.save(iris.cube.CubeList([s]),"/work/scratch-nopw/emmah/cmems/%04d%02d/%04d%02d%02d_relax_7daymean.nc"%(year1,year2%100,t0.year,t0.month,t0.day))
#  return
  mask = (t.data >0)
  ct = iris.Constraint(time = lambda t: t.point.replace(hour=0,minute=0) in times)
  if year1>2006:
    T_o = iris.load("/gws/nopw/j04/terramaris/emmah/sst_products/ostia/%04d1023_%04d0322_ostia_sst.nc"%(year1,year2),ct)[0]
    T_o.coord('latitude').guess_bounds()
    T_o.coord('longitude').guess_bounds()
    T_o = T_o.regrid(t,iris.analysis.AreaWeighted())
    T_o.convert_units("degreesC")
    T_o.rename("ostia")
    T_o = T_o.collapsed("time",iris.analysis.MEAN)
    T_o.coord('time').rename('t')
  else:
    try:
      T_o = iris.load("/gws/nopw/j04/terramaris/emmah/cmems_ostia/hadisst/era5_sst_%04d_?.nc"%year1,ct)
      equalise_attributes(T_o)
      T_o = T_o.concatenate_cube()
    except iris.exceptions.ConcatenateError:
      T_o = iris.load("/gws/nopw/j04/terramaris/emmah/cmems_ostia/hadisst/era5_sst_%04d_?.nc"%year1)
      equalise_attributes(T_o)
      T_o = T_o.concatenate_cube().extract(ct)
    T_o.coord('latitude').guess_bounds()
    T_o.coord('longitude').guess_bounds()
    T_o = T_o.collapsed("time",iris.analysis.MEAN)
    ind = nd.distance_transform_edt(T_o.data.mask,return_distances=False,return_indices=True)
    T_o.data = T_o.data[tuple(ind)]
    T_o = T_o.regrid(t,iris.analysis.Linear())
    T_o.convert_units("degreesC")
    T_o.coord('time').rename('t')
  delta = T_o - t[0]
#  inv = delta.data
#  ind = nd.distance_transform_edt(inv.mask,return_distances=False,return_indices=True)
#  delta.data = inv[tuple(ind)]
  delta.data=delta.data.filled(0)
  delta.add_aux_coord(t[0].coord('z'))
  delta.add_aux_coord(t[0].coord('t'))
  new = t+delta
  new.data = new.data*mask
  new.rename("temperature")
  iris.save(iris.cube.CubeList([new,s]),"/work/scratch-nopw/emmah/cmems/%04d%02d/%04d%02d%02d_relax_7daymean.nc"%(year1,year2%100,t0.year,t0.month,t0.day))

def combine(i,year1):
  year2 = year1+1
  t0 = dt.datetime(year1,11,1)+dt.timedelta(6*i)
  times = [t0+dt.timedelta(i) for i in range(6+1)]
  os.chdir("/work/scratch-nopw/emmah/")
  data= iris.load(["/work/scratch-nopw/emmah/cmems/%04d%02d/%04d%02d%02d_relax_7daymean.nc"%(year1,year2%100,t.year,t.month,t.day) for t in times])
  data = data.merge()
  print(data)
  os.chdir("/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/")
  ts = (t0-dt.datetime(year1,1,1)).days
  iris.save(data.extract("temperature")[0],"/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/%d/%04d_%05d_temperature_relax.nc"%(year1,year1,ts),zlib=True)
  iris.save(data.extract("salinity")[0],   "/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/%d/%04d_%05d_salinity_relax.nc"%(year1,year1,ts),zlib=True)


cmd = sys.argv[1]
year = int(sys.argv[2])

try:
  job = int(os.environ["SLURM_ARRAY_TASK_ID"])-1
except:
  job = int(sys.argv[3])

if cmd == "prep":
  for i in range(5):
    prep(job*5+i,year)
elif cmd =="roll":
  for i in range(4,5):
    roll(job*5+i,year)
elif cmd == "combine":
  combine(job,year)


#for job in range(1,140):
#  print(job)
#  roll(job)

#combine()

