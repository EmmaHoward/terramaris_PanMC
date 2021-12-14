#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH -o /home/users/emmah/log/corr-%a.o
#SBATCH -e /home/users/emmah/log/corr-%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=48000
import matplotlib.pyplot as plt
import numpy as np
import iris
import iris.plot as iplt
import cmocean
from panMC import panMC
import iris.cube
from matplotlib.ticker import MaxNLocator
cy =iris.Constraint(latitude=lambda y: -20<=y<=20)
cx =iris.Constraint(longitude=lambda y: 85<=y<=160)
cz=iris.Constraint(pressure=lambda z: z>=100)
era5path="/gws/nopw/j04/terramaris/emmah/era5/"

def plot(years,scratchpath,figname):
  MC12=iris.load([scratchpath+"MC12_pressurelevel_%04d.nc"%year for year in years],cx&cy&cz)
  era5 = iris.load(["/gws/nopw/j04/terramaris/emmah/era5/*vert_%04d12.nc"%year for year in years]+\
                 ["/gws/nopw/j04/terramaris/emmah/era5/*vert_%04d0?.nc"%(year+1) for year in years],cx&cy)
  iris.util.equalise_attributes(MC12)
  iris.util.equalise_attributes(era5)
  era5=era5.concatenate()
  for cube in MC12:
    cube.coord("time").convert_units("days since 2003-01-01")
  MC12=MC12.merge()
  era5 = iris.cube.CubeList([cube.collapsed(["longitude"],iris.analysis.MEAN) for cube in era5])
  era5 = iris.cube.CubeList([cube.rolling_window("time",iris.analysis.MEAN,3)[::3] for cube in era5])
  var = [ 'air_temperature','lagrangian_tendency_of_air_pressure', 'specific_humidity', 'eastward_wind', 'northward_wind']
  vmax= [1,0.03,0.0009,3.5,1.8 ]
  fig,axs=plt.subplots(2,5,figsize=(12,5))
  fig.subplots_adjust(top=0.88,left=0.09,right=0.97,hspace=0.2,wspace=0.4)
  for i,v in enumerate(var):
    ticks =MaxNLocator(10)._raw_ticks(0,vmax[i])
    ax1=axs[0,i]#plt.subplot(2,5,i+1)
    ax1.set_title(v.replace("_"," "))
    if i==1:
      ax1.set_title("omega")
    a=iplt.contourf(MC12.extract(v)[0].collapsed("time",iris.analysis.STD_DEV),ticks,cmap=cmocean.cm.dense,axes=ax1)
    CS=iplt.contour(MC12.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax1)
    ax1.set_ylim(1000,100)
    ax1.clabel(CS, CS.levels, inline=True, fontsize=6)
    ax2=axs[1,i]#plt.subplot(2,5,i+6)
    a=iplt.contourf(era5.extract(v)[0].collapsed("time",iris.analysis.STD_DEV),ticks,cmap=cmocean.cm.dense,axes=ax2)
    CS=iplt.contour(era5.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax2)
    ax2.clabel(CS, CS.levels, inline=True, fontsize=6)
    fig.colorbar(a,ax=[ax1,ax2],orientation="horizontal",ticks=a.levels[::3])
    ax2.set_ylim(1000,100)
    if i==0:
      ax1.set_ylabel("MC12 \n Pressure (hPa)")
      ax2.set_ylabel("ERA-5 \n Pressure (hPa)")
    ax2.set_xlabel("Latitude")
  plt.savefig(figname)
  plt.show()
