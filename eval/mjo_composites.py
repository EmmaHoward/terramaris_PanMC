#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
#SBATCH -p short-serial
#SBATCH --array=[1-18]
#SBATCH -o /home/users/emmah/log/mjo%a.o
#SBATCH -e /home/users/emmah/log/mjo%a.e 
#SBATCH -t 01:00:00
#SBATCH --mem=48000
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/precip/"
MC2_path  = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"

import matplotlib.pyplot as plt
from iris.util import equalise_attributes
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




def calc(mjo_phase,years,scratchpath,MC,rmm=True):
  if rmm:
    mjo = read_rmm()
    phase = mjo["mjo phase"]
  else:
    mjo=pd.read_csv('/home/users/emmah/reading/ivar/mjo.csv',index_col=0,names=['x0','x1','mag'],parse_dates=True,dayfirst=True)
    phase = np.ceil(((np.arctan2(mjo.x0,-mjo.x1))/(np.pi/4))%8)*(mjo.mag>=1)
  data= iris.cube.CubeList()
  for year in years:
    print(year)
    dates = phase[dt.datetime(year,12,1):dt.datetime(year+1,2,28)]
    dates = np.array(dates[dates==mjo_phase].index)
    print(len(dates))
#    dates= [dt.datetime(year,12,1)]
    if len(dates)>0:
      dates =pd.to_datetime(dates)
      if MC == "MC12":
        data += panMC(year,MC,"rainfall").load_iris(dates,variables=["convective_rainfall_amount","stratiform_rainfall_amount"])
      elif MC =="MC2":
        ct = iris.Constraint(time = lambda t: dt.datetime(t.point.year,t.point.month,t.point.day) in dates)
        data += iris.load([MC2_path+"%04d%02d_hourly_coarsened_rainfall.nc"%(yr,mon) for (yr,mon) in [(year,12),(year+1,1),(year+1,2)]],ct)
      print(data[-1])
  for cube in data:
    cube.coord("time").convert_units(data[0].coord("time").units)
  iris.util.equalise_attributes(data)
  import pdb;pdb.set_trace()
  data = data.concatenate()
  data= iris.cube.CubeList([cube.collapsed("time",iris.analysis.MEAN) for cube in data])
  if rmm:
    iris.save(data,scratchpath+"%s_mjo_rmm_rain_phase%d.nc"%(MC,mjo_phase),zlib=True)
  else:
    iris.save(data,scratchpath+"%s_mjo_rain_phase%d.nc"%(MC,mjo_phase),zlib=True)

def load(year,MC,path):
  if MC == "MC12":
    data=iris.load([path+"%04d%02d_diurnal_precip.nc"%(year+int(month<6),month) for month in [12,1,2]])
  else:
    data = iris.load([MC2_path+"%04d%02d_hourly_coarsened_rainfall.nc"%(yr,mon) for (yr,mon) in [(year,12),(year+1,1),(year+1,2)]])
  equalise_attributes(data)
  data=data.concatenate()
  if MC=="MC12":
    P = (data.extract("convective_rainfall_amount")[0]+data.extract("stratiform_rainfall_amount")[0])*4*24
    #P = (data[0]+data[1]+data[2]+data[3])*4*24
    P.units="mm/day"
  else:
    P = (data[0])*24
  return P.collapsed("time",iris.analysis.MEAN)

def abs_mean(years,scratchpath,MC):
  P=iris.cube.CubeList()
  for year in years:
    print(year)
    P.append(load(year,MC,MC12_path))
    P[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  equalise_attributes(P)
  mean = P.merge_cube().collapsed("time",iris.analysis.MEAN) 
  iris.save(mean,scratchpath+"%s_mean_rain.nc"%MC)

def plot(scratchpath,figname,rmm=True):
  fig=plt.figure(figsize=(9,7))
  cmap=cmap_discretise("BrBG",11)
  mean2 = iris.load(scratchpath+"MC2_mean_rain.nc")[0]*4
  mean12 = iris.load(scratchpath+"MC12_mean_rain.nc")[0]
  axs = []
  for phase in range(1,9):
    if rmm:
      data2=iris.load(scratchpath+"MC2_mjo_rmm_rain_phase%d.nc"%phase)
      data12=iris.load(scratchpath+"MC12_mjo_rmm_rain_phase%d.nc"%phase)
    else:
      data2=iris.load(scratchpath+"MC2_mjo_rain_phase%d.nc"%phase)
      data12=iris.load(scratchpath+"MC12_mjo_rain_phase%d.nc"%phase)
    data12 = (data12.extract("convective_rainfall_amount")[0]+data12.extract("stratiform_rainfall_amount")[0])*4*24
    data2 = (data2.extract("stratiform_rainfall_amount")[0])*24*4
    data12.units = "mm/day"
    print(data2.summary(shorten=True))
    print(data2.coord('time'))
   # data = (data[0]+data[1])*24
   # if phase not in [3,4]:
   #   data=data*4
    ax=plt.subplot(4,4,phase+4*(phase>4),projection=ccrs.PlateCarree())
    ax.coastlines()
    a=iplt.pcolormesh(data2-mean2,vmin=-11,vmax=11,cmap=cmap)
    plt.xlim(90,155)
    plt.ylim(-15,15)
    plt.title("MC2 MJO Phase %d"%phase)
    axs.append(ax)
    ax=plt.subplot(4,4,phase+4+4*(phase>4),projection=ccrs.PlateCarree())
    ax.coastlines()
    a=iplt.pcolormesh(data12-mean12,vmin=-11,vmax=11,cmap=cmap)
    plt.xlim(90,155)
    plt.ylim(-15,15)
    plt.title("MC12 MJO Phase %d"%phase)
    axs.append(ax)
  cax=fig.add_axes([0.3,0.1,0.4,0.05])
  fig.colorbar(a,cax=cax,orientation="horizontal",ticks=np.arange(-10,11,4))
#  fig.colorbar(a,ax=axs,ticks=np.arange(-10,11,4))
  fig.subplots_adjust(bottom=0.2,left=0.03,right=0.97)
  plt.savefig(figname)
  plt.show()


if __name__=="__main__":
  years=[2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]
  if 1:#len(os.environ["SLURM_ARRAY_TASK_ID"])>0:
      i = 9+3# int(os.environ["SLURM_ARRAY_TASK_ID"])-1
      mjo_phase,MC = i%9, ["MC2","MC12"][i//9]
      if mjo_phase > 0:
        calc(mjo_phase,years,"/home/users/emmah/eval/",MC,True)
      else:
        abs_mean(years,"/home/users/emmah/eval/",MC)
  else:
    plot("/home/users/emmah/eval/","OMM_12.png",True)


