#
# Make enso difference map for multi enso
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

cz = iris.Constraint(pressure=850)
cz2 = iris.Constraint(pressure_level=850)
cy =iris.Constraint(latitude=lambda y: -20<=y<=20)
cx =iris.Constraint(longitude=lambda y: 85<=y<=160)
ct=iris.Constraint(time=lambda t:t.point.month in [12,1,2])



def load(years):
  ymon = []
  for year in years:
    ymon +=[(year,12),(year+1,1),(year+1,2)]
  ct = iris.Constraint(time=lambda t: (t.point.year,t.point.month) in ymon)
  era5 = iris.load("/gws/nopw/j04/terramaris/emmah/era5/era5_monthly_means.nc",ct&cz2)
#  era5 = iris.cube.CubeList([cube.rolling_window("time",iris.analysis.MEAN,3)[::3] for cube in era5]).merge_cube()
  era5 = iris.cube.CubeList([cube.collapsed('time',iris.analysis.MEAN) for cube in era5])
  data = iris.cube.CubeList()
  domain="MC2"
  for i,year in enumerate(years):
    data+=iris.load("/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/monmean_pl/%s_monmean_pl_%d.nc"%(domain,year),cz).extract(["eastward_wind","northward_wind"])
  for cube in data:
    cube.coord("time").convert_units("days since 2003-01-01")
  iris.util.equalise_attributes(data)
  data2 = data.merge().concatenate()
  data2 = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data2])
  data = iris.cube.CubeList()
  domain="MC12"
  for i,year in enumerate(years):
    print(i)
    data+=iris.load("/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monmean_pl/%s_monmean_pl_%d.nc"%(domain,year),cz).extract(["eastward_wind","northward_wind"])
  for cube in data:
    cube.coord("time").convert_units("days since 2003-01-01")
  iris.util.equalise_attributes(data)
  data12 = data.merge().concatenate()
  data12 = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data12])
  return era5,data12,data2

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
  u = data.extract("eastward_wind")
  v = data.extract("northward_wind")
  ax.coastlines()
  a=iplt.contourf(P[0]-P[1],np.linspace(-11,11,12),cmap="bwr_r",extend='both')
  q=iplt.quiver((u[0]-u[1])[::k,::k],(v[0]-v[1])[::k,::k],scale=100,pivot="mid",color="0.2",headwidth=5)
  plt.quiverkey(q,0.9,-0.2,5,"5m/s")
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.colorbar(a,ticks=np.arange(-10,11,4))
  plt.xlim(85,160)
  plt.ylim(-20,20)
  ax.coastlines()
  plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
  plt.xlim(85,160)
  plt.ylim(-20,20)
  return a

def plot(years1,years2,figpath,fig):
  era51,MC121,MC21=load(years1)
  era52,MC122,MC22=load(years2)
  era5 = era51+era52
  MC2 = MC21+MC22
  MC12 = MC121+MC121
  P121,P21,Pref1 = precip_bias.load_all(years1,years1)
  P122,P22,Pref2 = precip_bias.load_all(years2,years2)
  import pdb;pdb.set_trace()
  P12=iris.cube.CubeList([P121.collapsed("time",iris.analysis.MEAN),P122.collapsed("time",iris.analysis.MEAN)]).merge_cube()
  P2=iris.cube.CubeList([P21.collapsed("time",iris.analysis.MEAN),P22.collapsed("time",iris.analysis.MEAN)]).merge_cube()
  Pref=iris.cube.CubeList([Pref1.collapsed("time",iris.analysis.MEAN),Pref2.collapsed("time",iris.analysis.MEAN)]).merge_cube()
#  plt.figure(figsize=(5,8))
  ax1=plt.subplot(322,projection=ccrs.PlateCarree())
  plt.title("(d) MC2")
  a=diff_plot(MC2,P2,ax1,20)
  ax2=plt.subplot(324,projection=ccrs.PlateCarree())
  plt.title("(e) MC12")
  a=diff_plot(MC12,P12,ax2,20)
  ax3=plt.subplot(326,projection=ccrs.PlateCarree())
  plt.title("(f) ERA5 + GPM-IMERG")
  a=diff_plot(era5,Pref,ax3,9)
  plt.savefig(figpath)
  plt.show()

fig=plt.figure()
plot([2007,2017],[2009,2015],'all_enso.png',fig)


"""
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
  plt.show()

"""
