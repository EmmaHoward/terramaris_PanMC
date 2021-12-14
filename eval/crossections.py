#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH -o /home/users/emmah/log/corr-%a.o
#SBATCH -e /home/users/emmah/log/corr-%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=48000
import matplotlib.pyplot as plt
import os
import numpy as np
import iris
import iris.plot as iplt
import cmocean
from panMC import panMC
import iris.cube
cy =iris.Constraint(latitude=lambda y: -20<=y<=20)
cx =iris.Constraint(longitude=lambda y: 85<=y<=160)

cz=iris.Constraint(pressure=lambda z: z>=100)
def load(scratchpath, years, MC, variables):
  xsections = iris.cube.CubeList()
  for year in years:
    out = iris.cube.CubeList()
    print(year)
    if 0:#os.path.exists("%s/%s_pressurelevel_%d.nc"%(scratchpath,MC,year)):
      out=iris.load("%s/%s_pressurelevel_%d.nc"%(scratchpath,MC,year),cz)
    else:
      H = panMC(year,MC,"Heavyside").load_iris()[0]
      for var in variables:# [ 'temperature_pd','omega', 'specific_humidity_pd', 'theta_w', 'eastward_wind_pd', 'northward_wind_pd']:
        data=panMC(year,MC,var).load_iris()[0]
        data.data.mask += H.data
        data = data.collapsed(["time","longitude"],iris.analysis.MEAN)
#      if MC == "MC2":
#        data=data.rolling_window("latitude",iris.analysis.MEAN,9)[::9]
        out.append(data)
        iris.save(out,"%s/%s_pressurelevel_%d_%s.nc"%(scratchpath,MC,year,var),zlib=True)
    xsections += out
  return xsections

def plot_single(data,ax,k):
  q = data.extract("specific_humidity")[0]
#  u = data.extract("eastward_wind")[0]
#  v = data.extract("northward_wind")[0]
  w = data.extract("lagrangian_tendency_of_air_pressure")[0]
#  plt.subplot(111,aspect=1/20)
  a=iplt.contourf(w,np.linspace(-0.11,0.11,12),cmap="BrBG_r")
  b=iplt.contour(q,np.arange(0,0.021,0.002),cmap=cmocean.cm.dense,zorder=2)
  c1=plt.colorbar(b)
  c1.ax.tick_params(labelsize=6)
  c2=plt.colorbar(a,ticks=np.arange(-0.1,0.11,0.04))
  c2.ax.tick_params(labelsize=6)
#  CS=iplt.contour(u,np.arange(-22,23,4),colors=["r"],linewidths=[1],zorder=3)
#  ax.clabel(CS, CS.levels, inline=True, fontsize=5)
  plt.ylim(1000,100)
  plt.grid(zorder=2)
  plt.xlabel("Longitude")
  plt.ylabel("Pressure (hPa)")
  lat = w.coord("latitude").points
#  try:
#    p = w.coord("pressure").points
#  except:
#    p = w.coord("pressure_level").points
#  q=plt.quiver(lat[::k],p,v[:,::k].data,-w[:,::k].data*20,pivot="mid",scale=70,headwidth=4,zorder=4,width=0.005)
  return q

def plot_diff(data,ax):
  q = data.extract("specific_humidity")[0]
#  u = data.extract("eastward_wind")[0]
#  v = data.extract("northward_wind")[0]
  w = data.extract("lagrangian_tendency_of_air_pressure")[0]
#  plt.subplot(111,aspect=1/20)
  a=iplt.contourf(w,np.arange(-0.045,0.046,0.01),cmap="BrBG_r",zorder=1)
  b=iplt.contour(q,np.arange(-0.00275,0.0028,0.0005),cmap="bwr_r",lws=[1],zorder=3)
  c1=plt.colorbar(b)
  c1.ax.tick_params(labelsize=6)
  c2=plt.colorbar(a,ticks=np.arange(-0.04,0.041,0.02))
  c2.ax.tick_params(labelsize=6)
  plt.ylim(1000,100)
  plt.grid(zorder=2)
 

def eq_diff(data,era5):
  new = iris.cube.CubeList()
  diff = iris.cube.CubeList()
  for cube in data:
    if len(era5.extract(cube.name()))>0:
       template = era5.extract(cube.name())[0]
       new.append(cube.interpolate([["latitude",template.coord("latitude").points]],iris.analysis.Linear()))
       diff.append(new[-1].copy(data = new[-1].data - template.data))
  return new,diff


def load_era5(years):
  out = iris.cube.CubeList()
  for year in years:
    out += iris.load(["/gws/nopw/j04/terramaris/emmah/era5/*vert_%04d%02d.nc"%(yr,mon) for (yr,mon) in [(year,12),(year+1,1),(year+1,2)]],cy)
  for cube in out:
    cube.coord("pressure_level").rename("pressure")
  return out

def plot(figname,MC12_years,MC2_years,scratchpath):
  MC12 = load(scratchpath,MC12_years,"MC12",['specific_humidity_pd','omega', 'eastward_wind_pd'])
  MC2 = load(scratchpath,MC2_years,"MC2",[ 'specific_humidity_pd', 'omega', 'eastward_wind_pd'])
  era5 = load_era5(MC2_years)
  iris.util.equalise_attributes(MC12)
  iris.util.equalise_attributes(MC2)
  iris.util.equalise_attributes(era5)
  era5=era5.concatenate()

  MC12 = MC12.extract(["lagrangian_tendency_of_air_pressure","specific_humidity"])
  MC2 = MC2.extract(["lagrangian_tendency_of_air_pressure","specific_humidity"])
  era5 = era5.extract(["lagrangian_tendency_of_air_pressure","specific_humidity"])

  print(era5)
  for cube in MC12+MC2:
    cube.coord("time").convert_units("days since 2003-01-01")

  MC12=MC12.merge()
  MC2=MC2.merge()

  era5 = iris.cube.CubeList([cube.collapsed(["time","longitude"],iris.analysis.MEAN) for cube in era5])
  if len(MC12_years)>1:
    MC12 = iris.cube.CubeList([cube.collapsed(["time"],iris.analysis.MEAN) for cube in MC12])
  if len(MC2_years)>1:
    MC2 = iris.cube.CubeList([cube.collapsed(["time"],iris.analysis.MEAN) for cube in MC2])

  for cube in MC2:
    cube = cube.rolling_window("latitude",iris.analysis.MEAN,9)

  MC12,diff_12=eq_diff(MC12,era5)
  MC2,diff_2=eq_diff(MC2,era5)
  diff=eq_diff(MC2,MC12)[1]

  fig=plt.figure(figsize=(10,7))
  ax=plt.subplot(331,aspect=1/20)
  ax.set_title("MC12")
  q1=plot_single(MC2,ax,8)
  ax=plt.subplot(335,aspect=1/20)
  ax.set_title("MC12")
  q2=plot_single(MC12,ax,8)
  ax=plt.subplot(339,aspect=1/20)
  ax.set_title("era-5")
  q3=plot_single(era5,ax,8)

  ax=plt.subplot(332,aspect=1/20)
  ax.set_title("MC2 - MC12")
  plot_diff(diff,ax)

  ax=plt.subplot(333,aspect=1/20)
  ax.set_title("MC2 - era-5")
  plot_diff(diff_2,ax)

  ax=plt.subplot(336,aspect=1/20)
  ax.set_title("MC12 - era-5")
  plot_diff(diff_12,ax)

  fig.tight_layout()
  plt.savefig(figname)
  plt.show()

