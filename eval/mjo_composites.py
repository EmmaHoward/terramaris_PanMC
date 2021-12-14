#!/apps/contrib/jaspy/miniconda_envs/jaspy3.7/m3-4.5.11/envs/jaspy3.7-m3-4.5.11-r20181219/bin/python
#SBATCH -p test
#SBATCH --array=[1-8]
#SBATCH -o /home/users/emmah/log/mjo%a.o
#SBATCH -e /home/users/emmah/log/mjo%a.e 
#SBATCH -t 01:00:00
#SBATCH --mem=48000
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monthly_diurnal_precip_amounts/"

import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import equalise_attributes
import os
import iris
from panMC import panMC
import datetime as dt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import iris.plot as iplt
from routines import cmap_discretise


def read_rmm():
  data = pd.read_csv("/home/users/emmah/RMM1RMM2.74toRealtime.txt",sep="\s+")
  date = pd.to_datetime(data['year']*10000+data['month']*100+data['day'],format="%Y%m%d")
  data.index=date
  data['mjo phase'] = data['phase']
  data['mag'] = data['amplitude']
  data["mjo phase"][data["mag"]<1]=0
  mjo = data[['RMM1',"RMM2","mjo phase","mag"]]
  return mjo




def calc(mjo_phase,years,scratchpath):
#  mjo=pd.read_csv('/home/users/emmah/reading/ivar/mjo.csv',index_col=0,names=['x0','x1','mag'],parse_dates=True,dayfirst=True)
  mjo = read_rmm()
  phase = mjo["mjo phase"]
  data= iris.cube.CubeList()
  for year in years:
    print(year)
    dates = phase[dt.datetime(year,12,1):dt.datetime(year+1,2,28)]
    dates = np.array(dates[dates==mjo_phase].index)
    print(len(dates))
#    dates= [dt.datetime(year,12,1)]
    if len(dates)>0:
      dates =pd.to_datetime(dates)
      data += panMC(year,"MC12","rainfall").load_iris(dates,variables=["convective_rainfall_amount","stratiform_rainfall_amount"])
      print(data[-1])
  for cube in data:
    cube.coord("time").convert_units(data[0].coord("time").units)
  data = data.concatenate()
  data= iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data])
  iris.save(data,scratchpath+"MC12_mjo_rmm_rain_phase%d.nc"%mjo_phase,zlib=True)

def load(year,MC,path):
  data=iris.load([path+"%04d%02d_diurnal_precip.nc"%(year+int(month<6),month) for month in [12,1,2]])
  equalise_attributes(data)
  data=data.concatenate()
  if MC=="MC12":
    P = (data.extract("convective_rainfall_amount")[0])*4*24
    #P = (data[0]+data[1]+data[2]+data[3])*4*24
  else:
    P = (data[0]+data[1])*4*24
  P.units="mm/day"
  return P.collapsed("time",iris.analysis.MEAN)

def abs_mean(years,scratchpath):
  P12=iris.cube.CubeList()
  for year in years:
    print(year)
    P12.append(load(year,"MC12",MC12_path))
    P12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  mean = P12.merge_cube().collapsed("time",iris.analysis.MEAN) 
  iris.save(mean,scratchpath+"MC12_mean_rain.nc")

def plot(scratchpath,figname):
  fig=plt.figure(figsize=(10,5))
  cmap=cmap_discretise("BrBG",11)
  mean = iris.load(scratchpath+"MC12_mean_rain.nc")[0]*4*24
  for phase in range(1,9):
    data=iris.load(scratchpath+"MC12_mjo_rmm_rain_phase%d.nc"%phase)
    data = (data.extract("convective_rainfall_amount")[0])*4*24
   # data = (data[0]+data[1])*24
   # if phase not in [3,4]:
   #   data=data*4
    ax=plt.subplot(2,4,phase,projection=ccrs.PlateCarree())
    ax.coastlines()
    a=iplt.pcolormesh(data-mean,vmin=-11,vmax=11,cmap=cmap)
    plt.xlim(90,155)
    plt.ylim(-15,15)
    plt.title("MJO Phase %d"%phase)
  cax=fig.add_axes([0.3,0.1,0.4,0.1])
  fig.colorbar(a,cax=cax,orientation="horizontal",ticks=np.arange(-10,11,4))
  fig.subplots_adjust(bottom=0.2,left=0.03,right=0.97)
  plt.savefig(figname)
  plt.show()


if __name__=="__main__":
  years=[2003,2007,2014,2015,2016,2017,2018]
  if 0:#len(os.environ["SLURM_ARRAY_TASK_ID"])>0:
    for mjo_phase in range(1,9):
      calc(mjo_phase,years,"/work/scratch-pw2/emmah/eval/")
  else:
    plot("/work/scratch-pw2/emmah/eval/","RMM.png")


