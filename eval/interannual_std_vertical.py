#
# plot standard deviation of interannual atmospheric structure (figure 11)
#


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

# use data loader from crossesctions plot code (figure 6) as same data is used
from crossections import load,load_era5

#domain spatial constraints
cy =iris.Constraint(latitude=lambda y: -20<=y<=20)
cx =iris.Constraint(longitude=lambda y: 90<=y<=155)
cz=iris.Constraint(pressure=lambda z: z>=100)

# file paths
era5path="/gws/nopw/j04/terramaris/emmah/era5/"
finalpathMC12 = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monmean_pl/"
finalpathMC2 = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/monmean_pl/"

def plot(MC12_years,MC2_years,scratchpath,figname):
  # read data
  MC12 = load(finalpathMC12,MC12_years,"MC12")
  MC2 = load(finalpathMC2,MC2_years,"MC2")
  era5 = load_era5(MC12_years)
  iris.util.equalise_attributes(MC12)
  iris.util.equalise_attributes(MC2)
  iris.util.equalise_attributes(era5)
  # annual means
  MC2 = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in MC2])
  MC12 = iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in MC12])
  for cube in MC12+MC2:
    cube.coord("time").convert_units(MC12[0].coord("time").units)
  MC12 = MC12.merge()
  MC2 = MC2.merge()
  era5=era5.concatenate()
  # annual, zonal means
  era5 = iris.cube.CubeList([cube.rolling_window("time",iris.analysis.MEAN,3)[::3] for cube in era5])
  era5 = iris.cube.CubeList([cube.collapsed("longitude",iris.analysis.MEAN) for cube in era5])
  var = [ 'air_temperature','lagrangian_tendency_of_air_pressure', 'specific_humidity', 'eastward_wind', 'northward_wind']
  # plot limits (total and bias)
  vmax= [1,0.03,0.0009,3.5,1.8 ]
  vmax2= [0.11,0.01,0.0003,0.5,0.5 ]
  # make plots
  fig,axs=plt.subplots(3,5,figsize=(12,8),sharex=True,sharey=True,subplot_kw={'aspect':0.03})
  a = {}
  # loop over variables (columns in plot)
  for i,v in enumerate(var):
    # colour ticks
    ticks =MaxNLocator(10)._raw_ticks(0,vmax[i])
    ticks2 =MaxNLocator(10)._raw_ticks(-vmax2[i],vmax2[i])
    # compute standard deviations
    std_2 = MC2.extract(v)[0].collapsed("time",iris.analysis.STD_DEV)
    std_12 = MC12.extract(v)[0].collapsed("time",iris.analysis.STD_DEV)
    std_e = era5.extract(v)[0].collapsed("time",iris.analysis.STD_DEV)
    # regrid for difference calculation
    std_2 = std_2.copy(data=std_2.data - std_e.interpolate([('latitude',std_2.coord('latitude').points)],iris.analysis.Linear()).data)
    std_12 = std_12.copy(data=std_12.data - std_e.interpolate([('latitude',std_12.coord('latitude').points)],iris.analysis.Linear()).data)
    # plot MC2 biases
    ax1=axs[1,i]#plt.subplot(2,5,i+1)
    ax1.set_title("("+'fghij'[i]+") "+v.replace("_"," "))
    a[i]=iplt.contourf(std_2,ticks2,cmap='PiYG',axes=ax1)
    cbar=fig.colorbar(a[i],ax=ax1,orientation='horizontal',ticks=a[i].levels[::4])
    if i==1:
      ax1.set_title('abcdefghijklmnop'[3*i+1]+"omega")
    CS=iplt.contour(MC2.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax1)
    ax1.set_ylim(1000,100)
    ax1.clabel(CS, CS.levels, inline=True, fontsize=6)
    # plot MC12 biases
    ax2=axs[2,i]#plt.subplot(2,5,i+1)
    ax2.set_title("("+'klmno'[i]+") "+v.replace("_"," "))
    a[i]=iplt.contourf(std_12,ticks2,cmap=f'PiYG',axes=ax2)
    #cbar=fig.colorbar(a[i],ax=ax2,orientation='horizontal',ticks=a[i].levels[::4])
    if i==1:
      ax2.set_title('abcdefghijklmnop'[3*i+1]+"omega")
    CS=iplt.contour(MC12.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax2)
    ax2.set_ylim(1000,100)
    ax2.clabel(CS, CS.levels, inline=True, fontsize=6)
    # plot era5 full fields
    ax3=axs[0,i]#plt.subplot(2,5,i+6)
    a[i]=iplt.contourf(era5.extract(v)[0].collapsed("time",iris.analysis.STD_DEV),ticks,cmap=cmocean.cm.dense,axes=ax3)
    fig.colorbar(a[i],ax=ax3,orientation='horizontal',ticks=a[i].levels[::4])
    CS=iplt.contour(era5.extract(v)[0].collapsed("time",iris.analysis.MEAN),colors="k",axes=ax3)
    ax3.clabel(CS, CS.levels, inline=True, fontsize=6)
    ax3.set_ylim(1000,100)
    ax3.set_title("("+'abcde'[i]+") "+v.replace("_"," "))
    if i==1:
      ax3.set_title("omega")
    if i==0:
      ax1.set_ylabel("MC12 - ERA5 \n Pressure (hPa)")
      ax2.set_ylabel("MC2 - ERA5 \n Pressure (hPa)")
      ax3.set_ylabel("ERA5 \n Pressure (hPa)")
    ax2.set_xlabel("Latitude")
    #cax=fig.add_axes([ax1.get_position().x0,0.03,ax1.get_position().width,0.05])
    #fig.colorbar(a[i],cax=cax,orientation="horizontal",ticks=a[i].levels[::3])
  plt.tight_layout()
  plt.savefig(figname)
