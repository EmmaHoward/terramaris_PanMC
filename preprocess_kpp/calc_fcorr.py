#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#SBATCH -p short-serial-4hr
#  SBATCH -p high-mem
#SBATCH -A short4hr
#SBATCH -o /home/users/emmah/log/fc_relax__.o
#SBATCH -e /home/users/emmah/log/fc_relax__.e 
#SBATCH -t 02:00:00
#SBATCH --mem=48000
import numpy as np
import os 
import iris
from iris.coord_categorisation import add_day_of_year
import datetime as dt
#job = int(os.environ["SLURM_ARRAY_TASK_ID"]) -1
#month =  [11,12,1,2,3][job]

def combine(month):
  years = [2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017]

  out=iris.cube.CubeList()
  for year in years:
    if month == 2 and year in [2003,2007,2011,2015]:
      ct = iris.Constraint(time=lambda t: t.point.day <=28)
      data=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/%04d/KPPocean_relaxrun_daymean_%04d%02d.nc"%(year,year+int(month<6),month),["heat_flux_correction_per_metre","virtual_salt_flux_correction_per_metre"]).extract(ct)
    elif month == 3 and year in [2003,2007,2011,2015]:
      ct = iris.Constraint(time=lambda t: t.point.day ==29)
      data1=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/%04d/KPPocean_relaxrun_daymean_%04d%02d.nc"%(year,year+int(month<6),month),["heat_flux_correction_per_metre","virtual_salt_flux_correction_per_metre"])
      data=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/%04d/KPPocean_relaxrun_daymean_%04d%02d.nc"%(year,year+int(month<6),2),["heat_flux_correction_per_metre","virtual_salt_flux_correction_per_metre"]).extract(ct)
      for cube in data1: 
        nt = cube.shape[0]
        for i in range(nt):
          data.append(cube[i])
      data = data.merge()  
    else:
      data=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/%04d/KPPocean_relaxrun_daymean_%04d%02d.nc"%(year,year+int(month<6),month),["heat_flux_correction_per_metre","virtual_salt_flux_correction_per_metre"])
    for cube in data:
      add_day_of_year(cube,"time","doyr")
      if year%4==0 and month > 6:
        cube.coord("doyr").points = cube.coord("doyr").points-1
      elif month < 6: 
        cube.coord("doyr").points = cube.coord("doyr").points+365
      cube.coord("time").convert_units("days since 2000-01-01")
      print(year,cube.coord("doyr").points)
      out.append(cube)
  out = out.concatenate()
  mean = iris.cube.CubeList()
  for cube in out:
    mean.append(cube.aggregated_by("doyr",iris.analysis.MEAN))
#    mean[-1].remove_coord("time")
  iris.save(mean,"/work/scratch-nopw/emmah/pp/fcorr_%02d.nc"%month)

def roll():
  ct = iris.Constraint(time=lambda t: t.point.day ==1)
  mask=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/2003/KPPocean_relaxrun_daymean_200311.nc","sea_water_temperature").extract(ct)[0].data==0
  data=iris.cube.CubeList()
  for month in [1,2,3]:
    tmp=iris.load("/work/scratch-nopw/emmah/pp/fcorr_%02d.nc"%month)
#    z=tmp[0].coord("depth").points
    data.append(tmp[0][:])#.copy(data=-1*np.trapz(tmp[0][:,:-1].data,z[:-1],axis=1)))
    data.append(tmp[1][:])#.copy(data=-1*np.trapz(tmp[1][:,:-1].data,z[:-1],axis=1)))
    t0 = (dt.datetime(2001+int(month<6),month,1)-dt.datetime(2000,12,31)).days
    print(t0)
    doyr = iris.coords.DimCoord(np.arange(t0+0.5,t0+{11:30,12:31,1:31,2:28,3:15}[month],1))
    doyr.rename("t")
#    if month==2:
#      data[-1]=data[-1][:-1]
#      data[-2]=data[-2][:-1]
    data[-1].remove_coord("time")
    data[-2].remove_coord("time")
    data[-1].add_dim_coord(doyr,0)
    data[-2].add_dim_coord(doyr,0)
    try:
      del data[-1].metadata.attributes["Conventions"]
      del data[-2].metadata.attributes["Conventions"]
    except:
      1
    for i in [-1,-2]:
      try:
        data[i].remove_coord("doyr")
      except:
        1
  data=data.concatenate()
  print(data)
#  data[0]=data[0][:,0]#:3,:3]
#  data[1]=data[1][:,0]#:3,:3]
  roll = [data[i][:].rolling_window('t',iris.analysis.MEAN,30) for i in range(2)]
  roll[0].coord("t").bounds=None
  roll[1].coord("t").bounds=None
  roll[0].coord("t").points = roll[0].coord("t").points - 0.5
  roll[1].coord("t").points = roll[1].coord("t").points - 0.5
#  start = iris.cube.CubeList()
#  for i in range(305,320):
#    start.append(data[0][:i-304+15].collapsed("t",iris.analysis.MEAN))
#    start.append(data[1][:i-304+15].collapsed("t",iris.analysis.MEAN))
#    start[-1].coord("t").points=i-0.5
#    start[-2].coord("t").points=i-0.5
#    start[-1].coord("t").bounds=None
#    start[-2].coord("t").bounds=None
#  start = start.merge()
  out = roll#iris.cube.CubeList(roll)#+start)
  print(out)
#  for cube in out:
#    cube.data = cube.data.astype("float")
#  out=out.concatenate()
#  iris.save(out,"/work/scratch-nopw/emmah/pp/fcorr_vint.nc")
  for i in range(382,420,6):
     ct = iris.Constraint(t=lambda tt: i<=tt<=i+9)
     hc = out.extract("heat_flux_correction_per_metre").extract(ct)[0]
     sc = out.extract("virtual_salt_flux_correction_per_metre").extract(ct)[0]
     hc.data.mask += mask
     sc.data.mask += mask
     hc.data = hc.data.filled(0)
     sc.data = sc.data.filled(0)
     hc.rename("fcorr")
     sc.rename("sfcorr")
     hc.coord("depth").rename("z")
     sc.coord("depth").rename("z")
     hc.coord("t").units="days since 1900-01-01"
     sc.coord("t").units="days since 1900-01-01"
     iris.save(sc,"/gws/nopw/j04/terramaris/emmah/coupled_N1280/HS_corrections/%05d_tm_15yr_scorr.nc"%i,zlib=True)
     iris.save(hc,"/gws/nopw/j04/terramaris/emmah/coupled_N1280/HS_corrections/%05d_tm_15yr_hcorr.nc"%i,zlib=True)

roll()
#combine(month)
