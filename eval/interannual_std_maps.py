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


compute_seas_mean(2007,"MC12")

#for year in years[:1]:
#  compute_seas_mean(year,"MC12")


def load(MC12_years,MC2_years):
  ymon = []
  for year in MC12_years:
    ymon +=[(year,12),(year+1,1),(year+1,2)]
  ct = iris.Constraint(time=lambda t: (t.point.year,t.point.month) in ymon)
  era5 = iris.load("/gws/nopw/j04/terramaris/emmah/era5/era5_monthly_means.nc",ct&cz2)
  era5 = iris.cube.CubeList([cube.rolling_window("time",iris.analysis.MEAN,3)[::3] for cube in era5])
  data = iris.cube.CubeList()
  domain="MC2"
  for i,year in enumerate(MC2_years):
    data+=iris.load("/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/monmean_pl/%s_monmean_pl_%d.nc"%(domain,year),cz).extract(["eastward_wind","northward_wind"])
  for cube in data:
    cube.coord("time").convert_units("days since 2003-01-01")
  data = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data])
  iris.util.equalise_attributes(data)
  data2 = data.merge()
  data = iris.cube.CubeList()
  domain="MC12"
  for i,year in enumerate(MC12_years):
    data+=iris.load("/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monmean_pl/%s_monmean_pl_%d.nc"%(domain,year),cz).extract(["eastward_wind","northward_wind"])
  for cube in data:
    cube.coord("time").convert_units("days since 2003-01-01")
  data = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data])
  iris.util.equalise_attributes(data)
  data12 = data.merge()
  return era5,data12,data2

def std_plot_1d(data,ax,std_max,mean_ticks,k=None):
  std = data.collapsed("time",iris.analysis.STD_DEV)
  mean = data.collapsed("time",iris.analysis.MEAN)
  if not k is None:
    if mean.data.mask.sum() > 10000:
       tmp = mean.data
       tmp[54:-55,42:-42] = gaussian_filter(mean[54:-55,42:-42].data,sigma=k)
       mean.data = np.ma.masked_values(tmp,0)
    else:
      mean.data = gaussian_filter(mean.data,sigma=k)
#  try:
#    std.data=std.data.filled(0)
#  except:
#    1
  ax.coastlines()
  a=iplt.contourf(std,np.linspace(0,std_max,11),cmap=cmocean.cm.dense,extend='max')
  CS = iplt.contour(mean,mean_ticks,colors="0.2",linewidths=1)
  ax.clabel(CS, CS.levels, inline=True, fontsize=5)
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a


def std_plot_2d(data,ax,k):
  u = data.extract("eastward_wind")[0]
  v = data.extract("northward_wind")[0]
  std = (u.collapsed("time",iris.analysis.VARIANCE)+v.collapsed("time",iris.analysis.VARIANCE))**0.5
  um = u.collapsed("time",iris.analysis.MEAN)
  vm = v.collapsed("time",iris.analysis.MEAN)
  try:
    std.data=std.data.filled(0)
  except:
    1
  ax.coastlines()
  a=iplt.contourf(std,np.linspace(0,3.5,15),cmap=cmocean.cm.dense,extend='max')
  q=iplt.quiver(um[::k,::k],vm[::k,::k],scale=200,pivot="mid",color="0.2")
  plt.quiverkey(q,0.9,-0.2,5,"5m/s")
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a

def main(MC12_years,MC2_years,figname):
  P12,P2,Pref = precip_bias.load_all(MC12_years,MC2_years)

  add_year(P2,"time","year")
  add_month(P2,"time","month")
  P2 = P2.aggregated_by(["month","year"],iris.analysis.MEAN)
  P2 = P2.rolling_window("time",iris.analysis.MEAN,3)[::3]
 
  add_year(P12,"time","year")
  add_month(P12,"time","month")
  P12 = P12.aggregated_by(["month","year"],iris.analysis.MEAN)
  P12 = P12.rolling_window("time",iris.analysis.MEAN,3)[::3]
  Pref = Pref.extract(ct).rolling_window("time",iris.analysis.MEAN,3)[::3]

  sst_12,sst_2,sstref = sst_bias.load_all(MC12_years,MC2_years)
  add_year(sst_2,"time","year")
  sst_2 = sst_2.aggregated_by(["month","year"],iris.analysis.MEAN)
  sst_2 = sst_2.rolling_window("time",iris.analysis.MEAN,3)[::3]
  add_year(sst_12,"time","year")
  add_month(sst_12,"time","month")
  sst_12 = sst_12.aggregated_by(["month","year"],iris.analysis.MEAN)
  sst_12 = sst_12.rolling_window("time",iris.analysis.MEAN,3)[::3]
  sstref = sstref.extract(ct).rolling_window("time",iris.analysis.MEAN,3)[::3]

  era5,MC12,MC2=load(MC12_years,MC2_years)
  fig=plt.figure(figsize=(12,9))

  ax1=plt.subplot(331,projection=ccrs.PlateCarree())
  a=std_plot_2d(MC2,ax1,20)
  ax1.set_title("MC2 Winds")

  ax2=plt.subplot(332,projection=ccrs.PlateCarree())
  a=std_plot_2d(MC12,ax2,20)
  ax2.set_title("MC12 Winds")

  ax3=plt.subplot(333,projection=ccrs.PlateCarree())
  a=std_plot_2d(era5,ax3,9)
  ax3.set_title("ERA-5 Winds")

  ax4=plt.subplot(334,projection=ccrs.PlateCarree())
  b=std_plot_1d(sst_2[:,3],ax4,1,np.arange(20,32,1))
  ax4.set_title("MC2 SST")

  ax5=plt.subplot(335,projection=ccrs.PlateCarree())
  b=std_plot_1d(sst_12[:,3],ax5,1,np.arange(20,32,1))
  ax5.set_title("MC12 SST")

  ax6=plt.subplot(336,projection=ccrs.PlateCarree())
  b=std_plot_1d(sstref[:,3],ax6,1,np.arange(20,32,1))
  ax6.set_title("Reference SST")

  ax7=plt.subplot(337,projection=ccrs.PlateCarree())
  c=std_plot_1d(P2[:],ax7,10,np.arange(0,25,2),20)
  ax7.set_title("MC2 Precip")

  ax8=plt.subplot(338,projection=ccrs.PlateCarree())
  c=std_plot_1d(P12[:],ax8,10,np.arange(0,25,2),20)
  ax8.set_title("MC12 Precip")

  ax9=plt.subplot(339,projection=ccrs.PlateCarree())
  c=std_plot_1d(Pref[:],ax9,10,np.arange(0,25,2),20)
  ax9.set_title("GPM Precip")
  plt.subplots_adjust(top=0.98,bottom=0.02,left=0.03,right=0.9,hspace=0.0,wspace=0.1)

  cax1 = fig.add_axes([0.92,ax1.get_position().y0,0.04,ax1.get_position().height])
  cax2 = fig.add_axes([0.92,ax4.get_position().y0,0.04,ax4.get_position().height])
  cax3 = fig.add_axes([0.92,ax7.get_position().y0,0.04,ax7.get_position().height])

  fig.colorbar(a,cax=cax1)#,orientation="horizontal")
  fig.colorbar(b,cax=cax2)#,orientation="horizontal")
  fig.colorbar(c,cax=cax3)#,orientation="horizontal")
  fig.savefig(figname)
  plt.show()


