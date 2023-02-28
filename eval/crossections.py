#
# Create cross-section plots (figure 6)
#


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
cy =iris.Constraint(latitude=lambda y: -15<=y<=15)
cx =iris.Constraint(longitude=lambda y: 90<=y<=155)
cz=iris.Constraint(pressure=lambda z: z>=100)

path= {"MC12": "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monmean_pl/",
       "MC2":  "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/monmean_pl/"}

def load(scratchpath, years, MC):
  # load monthly mean pressure level data
  xsections = iris.cube.CubeList()
  for year in years:
    print(year)
    data = iris.load(path[MC]+"%s_monmean_pl_%d.nc"%(MC,year),cx&cy&cz)
    for cube in data:# [ 'temperature_pd','omega', 'specific_humidity_pd', 'theta_w', 'eastward_wind_pd', 'northward_wind_pd']:
        cube.coord('time').convert_units('days since 2003-12-01T0000Z')
        cube = cube.collapsed("longitude",iris.analysis.MEAN)
        cube.data=cube.data.astype(np.float)
        xsections.append(cube)
  iris.util.equalise_attributes(xsections)
  xsections = xsections.merge()
  return xsections

def plot_single(data,ax,k):
  # plot mean state of data on axes ax
  q = data.extract("specific_humidity")[0]
  u = data.extract("eastward_wind")[0]
  w = data.extract("lagrangian_tendency_of_air_pressure")[0]
  a=iplt.contourf(w,np.linspace(-0.11,0.11,12),cmap="BrBG_r",extend='both')
  b=iplt.contour(q*1000,1000*np.arange(0,0.021,0.002),colors='k',linewidths=[1],zorder=2)
  c2=plt.colorbar(a,ticks=np.arange(-0.1,0.11,0.04))
  c2.ax.tick_params(labelsize=6)
  ax.clabel(b, b.levels, inline=True, fontsize=6,colors='k')
  plt.ylim(1000,100)
  plt.grid(zorder=2)
  plt.xlabel("Longitude")
  plt.ylabel("Pressure (hPa)")
  lat = w.coord("latitude").points
  return q

def plot_diff(data,ax):
  # plot biases of difference data
  q = data.extract("specific_humidity")[0]
#  u = data.extract("eastward_wind")[0]
#  v = data.extract("northward_wind")[0]
  w = data.extract("lagrangian_tendency_of_air_pressure")[0]
#  plt.subplot(111)
  a=iplt.contourf(w,np.arange(-0.045,0.046,0.01),cmap="BrBG_r",zorder=1)
  b1=iplt.contour(q*1000,1000*np.arange(0.0001,0.0012,0.0002),colors="b",linewidths=[1],zorder=3)
  b2=iplt.contour(q*1000,1000*np.arange(-0.0011,0.0,0.0002),colors="r",linewidths=[1],zorder=3,linestyles='solid')
  ax.clabel(b1, b1.levels, inline=True, fontsize=6)
  ax.clabel(b2, b2.levels, inline=True, fontsize=6)
  #c1.ax.tick_params(labelsize=6)
  c2=plt.colorbar(a,ticks=np.arange(-0.04,0.05,0.02))
  c2.ax.tick_params(labelsize=6)
  plt.ylim(1000,100)
  plt.grid(zorder=2)
 

def eq_diff(data,era5):
  # equalise data for comparison and then compute differences
  new = iris.cube.CubeList()
  diff = iris.cube.CubeList()
  for cube in data:
    if len(era5.extract(cube.name()))>0:
       template = era5.extract(cube.name())[0]
       new.append(cube.interpolate([["latitude",template.coord("latitude").points]],iris.analysis.Linear()))
       diff.append(new[-1].copy(data = new[-1].data - template.data))
  return new,diff


def load_era5(years):
  # load era5 monthly mean data
  ymon = []
  for year in years:
    ymon +=[(year,12),(year+1,1),(year+1,2)]
  ct = iris.Constraint(time=lambda t: (t.point.year,t.point.month) in ymon)
  out = iris.load("/gws/nopw/j04/terramaris/emmah/era5/era5_monthly_means.nc",ct&cx&cy)
  for cube in out:
    cube.coord("pressure_level").rename("pressure")
  return out

def plot(figname,MC12_years,MC2_years,scratchpath):
  # load data and generate figure
  MC2 = load(path,MC2_years,"MC2")
  MC12 = load(path,MC12_years,"MC12")
  era5 = load_era5(MC12_years)
  iris.util.equalise_attributes(MC12)
  iris.util.equalise_attributes(MC2)
  iris.util.equalise_attributes(era5)
  era5=era5.concatenate()
  print(era5)
  for cube in MC12+MC2:
    cube.coord("time").convert_units("days since 2003-01-01")
  MC12=MC12.merge()
  MC2=MC2.merge()
  # compute zonal /time means
  era5 = iris.cube.CubeList([cube.collapsed(["time","longitude"],iris.analysis.MEAN) for cube in era5])
  #if len(MC12_years)>1:
  MC12 = iris.cube.CubeList([cube.collapsed(["time"],iris.analysis.MEAN) for cube in MC12])
  #if len(MC2_years)>1:
  MC2 = iris.cube.CubeList([cube.collapsed(["time"],iris.analysis.MEAN) for cube in MC2])
  # smooth MC2 data to MC12 resolution
  for cube in MC2:
    cube = cube.rolling_window("latitude",iris.analysis.MEAN,9)
  # calculate differences
  MC12,diff_12=eq_diff(MC12,era5)
  MC2,diff_2=eq_diff(MC2,era5)
  diff=eq_diff(MC2,MC12)[1]
  # make figures
  fig=plt.figure(figsize=(9,7))
  ax=plt.subplot(331)
  ax.set_title("MC2")
  q1=plot_single(MC2,ax,8)
  ax=plt.subplot(335)
  ax.set_title("MC12")
  q2=plot_single(MC12,ax,8)
  ax=plt.subplot(339)
  ax.set_title("era-5")
  q3=plot_single(era5,ax,8)
  ax=plt.subplot(332)
  ax.set_title("MC2 - MC12")
  plot_diff(diff,ax)
  ax=plt.subplot(333)
  ax.set_title("MC2 - era-5")
  plot_diff(diff_2,ax)
  ax=plt.subplot(336)
  ax.set_title("MC12 - era-5")
  plot_diff(diff_12,ax)
  fig.tight_layout()
  plt.savefig(figname)
  plt.show()
