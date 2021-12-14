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

cz = iris.Constraint(pressure=850)
cz2 = iris.Constraint(pressure_level=850)
cy =iris.Constraint(latitude=lambda y: -20<=y<=20)
cx =iris.Constraint(longitude=lambda y: 85<=y<=160)
ct=iris.Constraint(time=lambda t:t.point.month in [12,1,2])



def compute_seas_mean(year,domain):
  u = panMC(year,domain, "eastward_wind_pd").load_iris(Constraints=cz)[0]
  v = panMC(year,domain,"northward_wind_pd").load_iris(Constraints=cz)[0]
  u = u.collapsed("time",iris.analysis.MEAN)
  v = v.collapsed("time",iris.analysis.MEAN)
  iris.save([u,v],"/work/scratch-pw2/emmah/%s_850wind_%d.nc"%(domain,year),zlib=True)


#for year in years[:1]:
#  compute_seas_mean(year,"MC12")
compute_seas_mean(2007,"MC12")


def load(years):
  era5 = iris.load(["/gws/nopw/j04/terramaris/emmah/era5/*vert_%04d12.nc"%year for year in years]+\
                 ["/gws/nopw/j04/terramaris/emmah/era5/*vert_%04d0?.nc"%(year+1) for year in years],cx&cy&cz2).extract(["eastward_wind","northward_wind"])
  iris.util.equalise_attributes(era5)
  era5=era5.concatenate()
  era5 = iris.cube.CubeList([cube.rolling_window("time",iris.analysis.MEAN,3)[::3] for cube in era5])
  data = iris.cube.CubeList()
  domain="MC12"
  for i,year in enumerate(years):
    data+=iris.load("/work/scratch-pw2/emmah/%s_850wind_%d.nc"%(domain,year))
  for cube in data:
    cube.remove_coord("forecast_reference_time")
    cube.coord("time").convert_units("days since 2003-01-01")
  data = data.merge()
  return era5,data

def std_plot_1d(data,ax,std_max,mean_ticks,k=None):
  std = data.collapsed("time",iris.analysis.STD_DEV)
  mean = data.collapsed("time",iris.analysis.MEAN)
  if not k is None:
    mean.data = gaussian_filter(mean.data,sigma=k)
  try:
    std.data=std.data.filled(0)
  except:
    1
  ax.coastlines()
  a=iplt.contourf(std,np.linspace(0,std_max,11),cmap="cubehelix_r")
  CS = iplt.contour(mean,mean_ticks,colors="0.2",linewidths=1)
  ax.clabel(CS, CS.levels, inline=True, fontsize=5)
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a


def diff_plot(data,P,ax,k):
  u = data.extract("eastward_wind")[0]
  v = data.extract("northward_wind")[0]
  ax.coastlines()
  a=iplt.contourf(P[0]-P[1],np.linspace(-11,11,12),cmap="bwr_r")
  q=iplt.quiver((u[0]-u[1])[::k,::k],(v[0]-v[1])[::k,::k],scale=100,pivot="mid",color="0.2")
  plt.quiverkey(q,0.9,-0.2,1,"1m/s")
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.colorbar(a)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a


if 1:
  ax=plt.subplot(311,projection=ccrs.PlateCarree())
  ax.coastlines()
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  ax2=plt.subplot(312,projection=ccrs.PlateCarree())
  a=diff_plot(MC12,P12,ax2,20)
  ax3=plt.subplot(313,projection=ccrs.PlateCarree())
  a=diff_plot(era5,Pref,ax3,9)

  import pdb;pdb.set_trace()
  fig=plt.figure(figsize=(9,9))

  ax1=plt.subplot(321,projection=ccrs.PlateCarree())
  ax1.set_title("MC12 Winds")

  ax2=plt.subplot(322,projection=ccrs.PlateCarree())
  a=std_plot_2d(era5,ax2,9)
  ax2.set_title("ERA-5 Winds")

  ax3=plt.subplot(323,projection=ccrs.PlateCarree())
  b=std_plot_1d(sst_12[:,3],ax3,1,np.arange(20,32,1))
  ax3.set_title("MC12 SST")

  ax4=plt.subplot(324,projection=ccrs.PlateCarree())
  b=std_plot_1d(sstref[:,3],ax4,1,np.arange(20,32,1))
  ax4.set_title("Reference SST")

  ax5=plt.subplot(325,projection=ccrs.PlateCarree())
  c=std_plot_1d(P12[:],ax5,10,np.arange(0,25,2),20)
  ax5.set_title("MC12 Precip")

  ax6=plt.subplot(326,projection=ccrs.PlateCarree())
  c=std_plot_1d(Pref[:],ax6,10,np.arange(0,25,2),20)
  ax6.set_title("GPM Precip")
  plt.subplots_adjust(top=0.98,bottom=0.02,left=0.03,right=0.9,hspace=0.0,wspace=0.1)

  cax1 = fig.add_axes([0.92,ax1.get_position().y0,0.04,ax1.get_position().height])
  cax2 = fig.add_axes([0.92,ax3.get_position().y0,0.04,ax3.get_position().height])
  cax3 = fig.add_axes([0.92,ax5.get_position().y0,0.04,ax5.get_position().height])

  fig.colorbar(a,cax=cax1)#,orientation="horizontal")
  fig.colorbar(b,cax=cax2)#,orientation="horizontal")
  fig.colorbar(c,cax=cax3)#,orientation="horizontal")
  fig.savefig(figname)
  plt.show()



