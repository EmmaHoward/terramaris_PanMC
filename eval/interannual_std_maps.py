#
# Generate maps of interannual standard deviations (Figure 10)
#

import iris
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import iris.cube
import cmocean
from iris.coord_categorisation import add_year
import iris.plot as iplt
import cartopy.crs as ccrs
from panMC import panMC
from iris.coord_categorisation import add_day_of_year,add_month
import numpy as np
import sst_bias
import precip_bias

# loading constraints
cz = iris.Constraint(pressure=850)
cz2 = iris.Constraint(pressure_level=850)
cy =iris.Constraint(latitude=lambda y: -20<=y<=20)
cx =iris.Constraint(longitude=lambda y: 85<=y<=160)
ct=iris.Constraint(time=lambda t:t.point.month in [12,1,2])


def load(MC12_years,MC2_years):
  # load input data
  ymon = []
  for year in MC12_years:
    ymon +=[(year,12),(year+1,1),(year+1,2)]
  ct = iris.Constraint(time=lambda t: (t.point.year,t.point.month) in ymon)
  #load era5 data
  era5 = iris.load("/gws/nopw/j04/terramaris/emmah/era5/era5_monthly_means.nc",ct&cz2)
  # construct annual means
  era5 = iris.cube.CubeList([cube.rolling_window("time",iris.analysis.MEAN,3)[::3] for cube in era5])
  data = iris.cube.CubeList()
  domain="MC2"
  # load MC2 data
  for i,year in enumerate(MC2_years):
    data+=iris.load("/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/monmean_pl/%s_monmean_pl_%d.nc"%(domain,year),cz).extract(["eastward_wind","northward_wind"])
  for cube in data:
    cube.coord("time").convert_units("days since 2003-01-01")
  # construct annual  means
  data = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data])
  iris.util.equalise_attributes(data)
  data2 = data.merge()
  data = iris.cube.CubeList()
  # lad MC12 data
  domain="MC12"
  for i,year in enumerate(MC12_years):
    data+=iris.load("/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monmean_pl/%s_monmean_pl_%d.nc"%(domain,year),cz).extract(["eastward_wind","northward_wind"])
  for cube in data:
    cube.coord("time").convert_units("days since 2003-01-01")
  # construct annual means
  data = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data])
  iris.util.equalise_attributes(data)
  data12 = data.merge()
  return era5,data12,data2

def std_plot_1d_diff(data1,data2,ax,std_max,mean_ticks,k=None):
  # make difference plots of scalar data
  # data1: first dataset
  # data2: second dataset
  # ax: plotting axes
  # std_max: maximum value of standard deviation to use in plotting colour axes
  # mean_ticks: contour levels for plotting mean state contour map
  # k: smoothing parameter
  std1 = data1.collapsed("time",iris.analysis.STD_DEV)
  std2 = data2.collapsed("time",iris.analysis.STD_DEV)
  std2.coord('longitude').coord_system = std1.coord('longitude').coord_system
  std2.coord('latitude').coord_system = std1.coord('latitude').coord_system
  std = std1 - std2.regrid(std1,iris.analysis.Linear())
  mean = data1.collapsed("time",iris.analysis.MEAN)
  # apply smoothing to mean 
  if not k is None:
    if mean.data.mask.sum() > 10000: # if data is masked, handle the mask before smoothing
       tmp = mean.data
       tmp[54:-55,42:-42] = gaussian_filter(mean[54:-55,42:-42].data,sigma=k)
       mean.data = np.ma.masked_values(tmp,0)
    else:
      mean.data = gaussian_filter(mean.data,sigma=k)
  # construct plot
  ax.coastlines()
  a=iplt.contourf(std,np.linspace(-std_max,std_max,12),cmap='PiYG',extend='max')
  CS = iplt.contour(mean,mean_ticks,colors="0.2",linewidths=1)
  ax.clabel(CS, CS.levels, inline=True, fontsize=5)
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a


def std_plot_1d(data,ax,std_max,mean_ticks,k=None):
  # make full field plots of 1-dimensional data
  # data: dataset
  # ax: plotting axes
  # std_max: maximum value of standard deviation to use in plotting colour axes
  # mean_ticks: contour levels for plotting mean state contour map
  # k: smoothing parameter
  std = data.collapsed("time",iris.analysis.STD_DEV)
  mean = data.collapsed("time",iris.analysis.MEAN)
  # apply smoothing to mean
  if not k is None:
    if mean.data.mask.sum() > 10000:# if data is masked, handle the mask before smoothing
       tmp = mean.data
       tmp[54:-55,42:-42] = gaussian_filter(mean[54:-55,42:-42].data,sigma=k)
       mean.data = np.ma.masked_values(tmp,0)
    else:
      mean.data = gaussian_filter(mean.data,sigma=k)
#  construct plot
  ax.coastlines()
  a=iplt.contourf(std,np.linspace(0,std_max,11),cmap=cmocean.cm.dense,extend='max')
  CS = iplt.contour(mean,mean_ticks,colors="0.2",linewidths=1)
  ax.clabel(CS, CS.levels, inline=True, fontsize=5)
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a

def std_plot_2d_diff(data1,data2,ax,k):
  # make difference plots of vector data (i.e. winds - hard coded)
  # data1: first dataset
  # data2: second dataset
  # ax: plotting axes
  # k: number of grid-cells to skip when plotting quiver maps (to make 
  u = data2.extract("eastward_wind")[0]
  v = data2.extract("northward_wind")[0]
  # compute standard deviations as root mean of the u and v component variances
  std2 = (u.collapsed("time",iris.analysis.VARIANCE)+v.collapsed("time",iris.analysis.VARIANCE))**0.5
  u = data1.extract("eastward_wind")[0]
  v = data1.extract("northward_wind")[0]
  std1 = (u.collapsed("time",iris.analysis.VARIANCE)+v.collapsed("time",iris.analysis.VARIANCE))**0.5
  std2.coord('longitude').coord_system = std1.coord('longitude').coord_system
  std2.coord('latitude').coord_system = std1.coord('latitude').coord_system
  # delta standard deviations between models
  std = std1 - std2.regrid(std1,iris.analysis.Linear())
  # time means for quiver plot
  um = u.collapsed("time",iris.analysis.MEAN)
  vm = v.collapsed("time",iris.analysis.MEAN)
  try: 
    std.data=std.data.filled(0)
  except:
    1
  # construct plots
  ax.coastlines()
  a=iplt.contourf(std,np.linspace(-1.1,1.1,12),cmap='PiYG')
  q=iplt.quiver(um[::k,::k],vm[::k,::k],scale=200,pivot="mid",color="0.2")
  plt.quiverkey(q,0.9,-0.2,5,"5m/s")
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a



def std_plot_2d(data,ax,k):
  # make full field plots of vector data (i.e. winds - hard coded)
  # data: dataset dataset
  # ax: plotting axes
  # k: number of grid-cells to skip when plotting quiver maps (to make 
  u = data.extract("eastward_wind")[0]
  v = data.extract("northward_wind")[0]
  std = (u.collapsed("time",iris.analysis.VARIANCE)+v.collapsed("time",iris.analysis.VARIANCE))**0.5
  um = u.collapsed("time",iris.analysis.MEAN)
  vm = v.collapsed("time",iris.analysis.MEAN)
  try:
    std.data=std.data.filled(0)
  except:
    1
  # make plot
  ax.coastlines()
  a=iplt.contourf(std,np.linspace(0,3.5,15),cmap=cmocean.cm.dense,extend='max')
  q=iplt.quiver(um[::k,::k],vm[::k,::k],scale=200,pivot="mid",color="0.2")
  plt.quiverkey(q,0.9,-0.2,5,"5m/s")
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a

def main(MC12_years,MC2_years,figname,cached=True):
  if not cached:
      # load precipitation data into memory
      P12,P2,Pref = precip_bias.load_all(MC12_years,MC2_years)
      add_year(P2,"time","year")
      add_month(P2,"time","month")
      # compute means for each season
      P2 = P2.aggregated_by(["month","year"],iris.analysis.MEAN)
      P2 = P2.rolling_window("time",iris.analysis.MEAN,3)[::3]
      add_year(P12,"time","year")
      add_month(P12,"time","month")
      P12 = P12.aggregated_by(["month","year"],iris.analysis.MEAN)
      P12 = P12.rolling_window("time",iris.analysis.MEAN,3)[::3]
      Pref = Pref.extract(ct).rolling_window("time",iris.analysis.MEAN,3)[::3]
      # load sea surface temperature data into memory
      sst_12,sst_2,sstref = sst_bias.load_all(MC12_years,MC2_years)
      add_year(sst_2,"time","year")
      # compute means for each season
      sst_2 = sst_2.aggregated_by(["month","year"],iris.analysis.MEAN)
      sst_2 = sst_2.rolling_window("time",iris.analysis.MEAN,3)[::3]
      add_year(sst_12,"time","year")
      add_month(sst_12,"time","month")
      sst_12 = sst_12.aggregated_by(["month","year"],iris.analysis.MEAN)
      sst_12 = sst_12.rolling_window("time",iris.analysis.MEAN,3)[::3]
      sstref = sstref.extract(ct).rolling_window("time",iris.analysis.MEAN,3)[::3]
      # load wind data into memory (returns means for each season)
      era5,MC12,MC2=load(MC12_years,MC2_years)    
      # cache data
      iris.save(MC2,"/home/users/emmah/eval/clim_MC2_uv.nc")
      iris.save(MC12,"/home/users/emmah/eval/clim_MC12_uv.nc")
      iris.save(era5,"/home/users/emmah/eval/clim_era5_uv.nc")
      iris.save(P2,"/home/users/emmah/eval/clim_MC2_P.nc")
      iris.save(P12,"/home/users/emmah/eval/clim_MC12_P.nc")
      iris.save(Pref,"/home/users/emmah/eval/clim_obs_P.nc")
      iris.save(sst_2,"/home/users/emmah/eval/clim_MC2_sst.nc")
      iris.save(sst_12,"/home/users/emmah/eval/clim_MC12_sst.nc")
      iris.save(sstref,"/home/users/emmah/eval/clim_obs_sst.nc")
  else:
      # load data from cache
      MC2 = iris.load("/home/users/emmah/eval/clim_MC2_uv.nc")
      MC12 = iris.load("/home/users/emmah/eval/clim_MC12_uv.nc")
      era5 = iris.load("/home/users/emmah/eval/clim_era5_uv.nc")
      P2 = iris.load_cube("/home/users/emmah/eval/clim_MC2_P.nc")
      P12 = iris.load_cube("/home/users/emmah/eval/clim_MC12_P.nc")
      Pref = iris.load_cube("/home/users/emmah/eval/clim_obs_P.nc")
      sst_2 = iris.load_cube("/home/users/emmah/eval/clim_MC2_sst.nc")
      sst_12 = iris.load_cube("/home/users/emmah/eval/clim_MC12_sst.nc")
      sstref = iris.load_cube("/home/users/emmah/eval/clim_obs_sst.nc")
  # generate plots
  fig=plt.figure(figsize=(12,7))
  ax1=plt.subplot(332,projection=ccrs.PlateCarree())
  a=std_plot_2d_diff(MC2,era5,ax1,20)
  ax1.set_title("(b) MC2 - Ref Winds")
  fig.colorbar(a,ax=ax1,ticks=np.arange(-1,1.1,0.4))

  ax2=plt.subplot(333,projection=ccrs.PlateCarree())
  a=std_plot_2d_diff(MC12,era5,ax2,20)
  ax2.set_title("(c) MC12 - Ref Winds")
  fig.colorbar(a,ax=ax2,ticks=np.arange(-1,1.1,0.4))

  ax3=plt.subplot(331,projection=ccrs.PlateCarree())
  a=std_plot_2d(era5,ax3,9)
  ax3.set_title("(a) ERA-5 Winds")
  fig.colorbar(a,ax=ax3)

  ax4=plt.subplot(335,projection=ccrs.PlateCarree())
  b=std_plot_1d_diff(sst_2[:,3],sstref[:,3],ax4,0.55,np.arange(20,32,1))
  ax4.set_title("(e) MC2 - Ref SST")
  fig.colorbar(b,ax=ax4,ticks=np.arange(-0.5,0.6,0.2))

  ax5=plt.subplot(336,projection=ccrs.PlateCarree())
  b=std_plot_1d_diff(sst_12[:,3],sstref[:,3],ax5,0.55,np.arange(20,32,1))
  ax5.set_title("(f) MC12 - Ref SST")
  fig.colorbar(b,ax=ax5,ticks=np.arange(-0.5,0.6,0.2))

  ax6=plt.subplot(334,projection=ccrs.PlateCarree())
  b=std_plot_1d(sstref[:,3],ax6,1,np.arange(20,32,1))
  ax6.set_title("(d) Reference SST")
  fig.colorbar(b,ax=ax6)

  ax7=plt.subplot(338,projection=ccrs.PlateCarree())
  c=std_plot_1d_diff(P2[:],Pref,ax7,2.75,np.arange(0,25,2),20)
  ax7.set_title("(h) MC2 - Ref Precip")
  fig.colorbar(c,ax=ax7,ticks=np.arange(-2.5,2.6,1))

  ax8=plt.subplot(339,projection=ccrs.PlateCarree())
  c=std_plot_1d_diff(P12[:],Pref,ax8,2.75,np.arange(0,25,2),20)
  ax8.set_title("(i) MC12 - Ref Precip")
  fig.colorbar(c,ax=ax8,ticks=np.arange(-2.5,2.6,1))

  ax9=plt.subplot(337,projection=ccrs.PlateCarree())
  c=std_plot_1d(Pref[:],ax9,10,np.arange(0,25,2),20)
  ax9.set_title("(g) GPM Precip")
  fig.colorbar(c,ax=ax9)
  plt.tight_layout()
  fig.savefig(figname)
  plt.show()


