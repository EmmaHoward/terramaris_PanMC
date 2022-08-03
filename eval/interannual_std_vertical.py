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
from crossections import load,load_era5

cy =iris.Constraint(latitude=lambda y: -20<=y<=20)
cx =iris.Constraint(longitude=lambda y: 90<=y<=155)
cz=iris.Constraint(pressure=lambda z: z>=100)
era5path="/gws/nopw/j04/terramaris/emmah/era5/"
finalpathMC12 = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monmean_pl/"
finalpathMC2 = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/monmean_pl/"

def plot(MC12_years,MC2_years,scratchpath,figname):
  MC12 = load(finalpathMC12,MC12_years,"MC12")
  MC2 = load(finalpathMC2,MC2_years,"MC2")
  era5 = load_era5(MC12_years)
  iris.util.equalise_attributes(MC12)
  iris.util.equalise_attributes(MC2)
  iris.util.equalise_attributes(era5)
  import pdb;pdb.set_trace()
  MC2 = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in MC2])
  MC12 = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in MC12])
  for cube in MC12+MC2:
    cube.coord("time").convert_units(MC12[0].coord("time").units)
  MC12 = MC12.merge()
  MC2 = MC2.merge()
  era5=era5.concatenate()
  era5 = iris.cube.CubeList([cube.rolling_window("time",iris.analysis.MEAN,3)[::3] for cube in era5])
  era5 = iris.cube.CubeList([cube.collapsed("longitude",iris.analysis.MEAN) for cube in era5])
  var = [ 'air_temperature','lagrangian_tendency_of_air_pressure', 'specific_humidity', 'eastward_wind', 'northward_wind']
  vmax= [1,0.03,0.0009,3.5,1.8 ]
  fig,axs=plt.subplots(3,5,figsize=(12,8))
  fig.subplots_adjust(top=0.956,bottom=0.15,left=0.08,right=0.988,hspace=0.301,wspace=0.355)
  a = {}
  for i,v in enumerate(var):
    ticks =MaxNLocator(10)._raw_ticks(0,vmax[i])
    ax1=axs[1,i]#plt.subplot(2,5,i+1)
    ax1.set_title(v.replace("_"," "))
    if i==1:
      ax1.set_title("omega")
    a[i]=iplt.contourf(MC2.extract(v)[0].collapsed("time",iris.analysis.STD_DEV),ticks,cmap=cmocean.cm.dense,axes=ax1)
    CS=iplt.contour(MC2.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax1)
    ax1.set_ylim(1000,100)
    ax1.clabel(CS, CS.levels, inline=True, fontsize=6)
    ax2=axs[0,i]#plt.subplot(2,5,i+1)
    ax2.set_title(v.replace("_"," "))
    if i==1:
      ax2.set_title("omega")
    a[i]=iplt.contourf(MC12.extract(v)[0].collapsed("time",iris.analysis.STD_DEV),ticks,cmap=cmocean.cm.dense,axes=ax2)
    CS=iplt.contour(MC12.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax2)
    ax2.set_ylim(1000,100)
    ax2.clabel(CS, CS.levels, inline=True, fontsize=6)
    ax3=axs[2,i]#plt.subplot(2,5,i+6)
    a[i]=iplt.contourf(era5.extract(v)[0].collapsed("time",iris.analysis.STD_DEV),ticks,cmap=cmocean.cm.dense,axes=ax3)
    CS=iplt.contour(era5.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax3)
    ax3.clabel(CS, CS.levels, inline=True, fontsize=6)
    ax3.set_ylim(1000,100)
    if i==0:
      ax1.set_ylabel("MC12 \n Pressure (hPa)")
      ax2.set_ylabel("MC2 \n Pressure (hPa)")
      ax3.set_ylabel("ERA-5 \n Pressure (hPa)")
    ax3.set_xlabel("Latitude")
    cax=fig.add_axes([ax1.get_position().x0,0.03,ax1.get_position().width,0.05])
    fig.colorbar(a[i],cax=cax,orientation="horizontal",ticks=a[i].levels[::3])
  plt.savefig(figname)
