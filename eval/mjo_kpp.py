#
# Make MJO-KPP longitude-phase plots (Figure 14)
#
import iris
from routines import pcoloraxis
import pickle
from scipy.interpolate import interp1d
from panMC import panMC
from matplotlib.colors import ListedColormap
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import subprocess
from iris.coord_categorisation import add_day_of_year
import pandas as pd
from scipy.ndimage import convolve
import seaborn as sns
import datetime as dt


cy=iris.Constraint(latitude=lambda y:-15<=y<=2)
# data paths
MC12_years = [2003,2007,2014,2015,2016,2017,2018]
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/sst/"
MC2_path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/sst/"
ref_path = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/"
sst12 = iris.cube.CubeList()
sstref = iris.cube.CubeList()

def read_rmm():
  # read RMM indices
  data = pd.read_csv("/home/users/emmah/RMM1RMM2.74toRealtime.txt",sep="\s+")
  date = pd.to_datetime(data['year']*10000+data['month']*100+data['day'],format="%Y%m%d")
  data.index=date
  data['mjo phase'] = data['phase']
  data['mag'] = data['amplitude']
  data["mjo phase"][data["mag"]<1]=0
  mjo = data[['RMM1',"RMM2","mjo phase","mag"]]
  return mjo


def load(year,MC,path,scratchpath=None):
  # load diurnal sst data
  outname = "%s_%04d%02d_sst_diurnal_range.nc"%(MC,year,(year+1)%100)
  print(year)
  if 1:# os.path.exists(path+outname):
    print("exists")
    dcmean = iris.load(path+outname)
  else:
    data=panMC(year,MC,"sea_water_temperature").load_iris()[0]
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    data.coord("time").bounds = np.array([np.round(data.coord("time").points)-1,np.round(data.coord("time").points)]).T
    data.coord("time").points = np.round(data.coord("time").points)-0.5
    add_day_of_year(data,"time","doyr")
    dcmax = data[:,0].aggregated_by("doyr",iris.analysis.MAX)
    dcmean = data[:,3].aggregated_by("doyr",iris.analysis.MEAN)
    dcmin = data[:,0].aggregated_by("doyr",iris.analysis.MIN)
    dc = dcmax - dcmin
    dc.rename("sst_diurnal_range")
    iris.save([dc,dcmean],scratchpath+outname,zlib=True)
    subprocess.check_call("cp %s%s %s%s"%(scratchpath,outname,path,outname),shell=True)
    subprocess.check_call("rm %s%s"%(scratchpath,outname),shell=True)
  return dcmean

#for year in MC12_years:
#  load(year,"MC12",MC12_path)

def detrend_seas(data,n):
  # remove seasonal cycle from n years of linear data
  nt,nx = data.shape
  tmp = data.data.reshape(n,nt//n,nx)
  mean = convolve(tmp.mean(axis=0),np.ones((10,1))/10)
  tmp = tmp - mean
  return data.copy(data=tmp.reshape(data.shape))


def load_kpp(years,domain,path,scratchpath=None):
  # load kpp mixed layer depth data and SST data
  # this will not be possible from preliminary available data for review, however dataset is saved
  # to pickle file later in code, which will be accessable
  # years: list of years to load
  # domain: MC2 or MC12
  # path:location of postprocessed data
  hmix,sst_m,sst_r = iris.cube.CubeList(),iris.cube.CubeList(),iris.cube.CubeList()
  for year in years:
    # load mixed layer depth
    data = panMC(year,domain,"mixedlayerdepth").load_iris(Constraints=cy)[0]
    data.coord("time").bounds = np.array([np.round(data.coord("time").points)-1,np.round(data.coord("time").points)]).T
    data.coord("time").points = np.round(data.coord("time").points)-0.5
    add_day_of_year(data,"time","doyr")
    hmix.append(data[:].aggregated_by("doyr",iris.analysis.MAX).collapsed("latitude",iris.analysis.MEAN))
    if year==2015:
        hmix[-1].coord("time").points = hmix[-2].coord("time").points
        hmix[-1].coord("time").bounds = hmix[-2].coord("time").bounds
    # load sea surface temperature: daily means and diurnal ranges
    data= load(year,domain,path,scratchpath=scratchpath)
    sst_m.append(data.extract_cube("sea_water_temperature").extract(cy).collapsed("latitude",iris.analysis.MEAN))
    sst_r.append(data.extract_cube("sst_diurnal_range").extract(cy).collapsed("latitude",iris.analysis.MEAN))
    if year==2015:
        sst_m[-1].coord("time").points = sst_m[-2].coord("time").points
        sst_r[-1].coord("time").points = sst_r[-2].coord("time").points
        sst_m[-1].coord("time").bounds = sst_m[-2].coord("time").bounds
        sst_r[-1].coord("time").bounds = sst_r[-2].coord("time").bounds
  # merge datasets into single cube
  for cube in hmix + sst_m + sst_r:
     cube.coord("time").convert_units("days since 2003-01-01")
     cube.remove_coord("doyr")
  hmix = hmix.concatenate_cube()
  sst_m = sst_m.concatenate_cube()
  sst_r = sst_r.concatenate_cube()
  # detrend data
  sst_m = detrend_seas(sst_m,len(years))
  hmix = hmix -  hmix.collapsed("time",iris.analysis.MEAN)
  sst_r = sst_r - sst_r.collapsed("time",iris.analysis.MEAN)
  return hmix,sst_m,sst_r


def mjo_composites(phase,hmix,sst_m,sst_r,years):
  """
   average mixed layer depth, daily mean sst and diurnal sst range across MJO phases
   phase: MJO timeseries
   hmix: mixed layer depth
   sst_m: daily mean sst
   sst_r: sst diurnal range
   years: seasons to consider
  """
  dates = pd.concat([phase[dt.datetime(year,12,1):dt.datetime(year+1,2,28)] for year in years])
  nt,nx = hmix.shape
  hmix_mjo,sst_r_mjo, sst_m_mjo = np.zeros((8,nx)),np.zeros((8,nx)),np.zeros((8,nx))
  for mjo_phase in range(1,9):
    dates2 =dates[dates==mjo_phase].index.date
    ct = iris.Constraint(time=lambda t: dt.date(t.point.year,t.point.month,t.point.day) in dates2)
    hmix_mjo[mjo_phase-1] = hmix.extract(ct).collapsed("time",iris.analysis.MEAN).data
    sst_m_mjo[mjo_phase-1] = sst_m.extract(ct).collapsed("time",iris.analysis.MEAN).data
    sst_r_mjo[mjo_phase-1] = sst_r.extract(ct).collapsed("time",iris.analysis.MEAN).data
  return hmix_mjo,sst_m_mjo,sst_r_mjo

def read_precip(scratchpath,MC):
  # load precipitation anomalies for each MJO phase
  # overall mean rainfall. 
  mean = iris.load(scratchpath+MC+"_mean_rain.nc",cy)[0]*{"MC2":4,"MC12":1}[MC] # factors are required to fix units 
  mask = iris.load("/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/invariants/land_binary_mask.nc",cy)[0]
  P = []
  for phase in range(1,9):
    # load rainfall data for each MJO phase
    data=iris.load(scratchpath+MC+"_mjo_rmm_rain_phase%d.nc"%phase,cy)
    if MC=="MC12":
      data = (data.extract('convective_rainfall_amount')[0]+data.extract('stratiform_rainfall_amount')[0])*24*4
    if MC=="MC2":
      data = (data[0])*24*4
      mask = mask.regrid(data,iris.analysis.Nearest())
    data.units = mean.units
    cube = data-mean
    # compute anomaly
    P.append(np.ma.masked_array(data.data-mean.data,0*mask.data).mean(axis=0))
  return np.ma.masked_array(P)

def main(plottype,calc,scratchpath,MC12_years,MC2_years,figname=None):
  # load mjo index
  mjo = read_rmm()
  phase = mjo["mjo phase"]#np.ceil(((np.arctan2(mjo.x0,-mjo.x1))/(np.pi/4))%8)*(mjo.mag>=1)
  # load data for first year (in order to access data grids)
  hmix12,sst_m12,sst_r12 = load_kpp(MC12_years[:1],"MC12",MC12_path,scratchpath=scratchpath)
  hmix2,sst_m2,sst_r2 = load_kpp(MC12_years[:1],"MC2",MC2_path,scratchpath=scratchpath)
  if calc:
    # compute MJO composites for ocean parameters
    hmix_mjo12,sst_m_mjo12,sst_r_mjo12 = mjo_composites(phase,hmix12,sst_m12,sst_r12,MC12_years,scratchpath)
    pickle.dump((hmix_mjo12,sst_m_mjo12,sst_r_mjo12),open(scratchpath+"kpp_mjo_MC12.pickle","wb"))
    hmix_mjo2,sst_m_mjo2,sst_r_mjo2 = mjo_composites(phase,hmix2,sst_m2,sst_r2,MC2_years,scratchpath)
    pickle.dump((hmix_mjo2,sst_m_mjo2,sst_r_mjo2),open(scratchpath+"kpp_mjo_MC2.pickle","wb"))
  else:
    # used cached MJO composites for ocean parameters
    (hmix_mjo12,sst_m_mjo12,sst_r_mjo12)=pickle.load(open(scratchpath+"kpp_mjo_MC12.pickle","rb"))
    (hmix_mjo2,sst_m_mjo2,sst_r_mjo2)=pickle.load(open(scratchpath+"kpp_mjo_MC2.pickle","rb"))
  # load rainfall MJO composites
  P_mjo12=read_precip(scratchpath,'MC12')
  P_mjo2=read_precip(scratchpath,"MC2")
  if not figname:
    return
  def smooth(data,n=50):
    # smoothing function
    return convolve(data,np.ones((1,n))/n)
  # start figure
  fig=plt.figure(figsize=(10,7))
  if plottype=="lines":
    # plot style 1: line plot (Not used)
    lon = hmix12.coord("longitude").points
    with sns.color_palette("turbo", 8):
      plt.subplot(141)
      plt.plot(lon,smooth(sst_r_mjo12[:]).T,lw=1)
      plt.xlabel("Longitude")
      plt.title("MC12 Diurnal SST Range")
      plt.subplot(142)
      a=plt.plot(lon,smooth(sst_m_mjo12[:]).T,lw=1)
      plt.title("MC12 Foundational SST")
      plt.xlabel("Longitude")
      plt.subplot(143)
      plt.plot(lon[10:-10],smooth(P_mjo12[:,10:-10]).T,lw=1)
      plt.xlabel("Longitude")
      plt.title("MC12 Precip over Ocean")
      plt.subplot(144)
      plt.plot(lon,smooth(hmix_mjo12[:]).T,lw=1)
      plt.xlabel("Longitude")
      plt.title("MC12 Mixed Layer Depth")
      fig.legend(a,["phase %d"%i for i in range(1,9,1)],loc="lower center",ncol=4)
    fig.subplots_adjust(top=0.85,bottom=0.3,left=0.06,right=0.97,hspace=0.2,wspace=0.246)

  elif plottype=="contour":
    # plot style 2: filled contour plot (Not used)
    indices = [-1,0,1,2,3,4,5,6,7,0]
    lon = hmix12.coord("longitude").points
    plt.subplot(141)
    plt.ylabel("MJO Phase")
    plt.grid()
    plt.contourf(lon,np.arange(10),smooth(sst_r_mjo12)[indices],np.linspace(-0.225,0.225,10),cmap="RdGy_r",extend="both")
    plt.xlabel("Longitude")
    plt.title("MC12 Diurnal SST Range")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-0.2,0.21,0.1))
    plt.subplot(142)
    a=plt.contourf(lon,np.arange(10),smooth(sst_m_mjo12)[indices],np.linspace(-0.45,0.45,10),cmap="bwr",extend="both")
    plt.grid()
    plt.title("MC12 Foundational SST")
    plt.xlabel("Longitude")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-0.4,0.5,0.2))
    plt.subplot(143)
    plt.contourf(lon[10:-10],np.arange(10),smooth(P_mjo12[:,10:-10])[indices],np.linspace(-9,9,10),cmap="BrBG",extend="both")
    plt.grid()
    plt.xlabel("Longitude")
    plt.title("MC12 Precip over Ocean")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-8,9,2))
    plt.subplot(144)
    plt.contourf(lon,np.arange(10),smooth(hmix_mjo12)[indices],np.linspace(-11,11,12),cmap="PiYG",extend="both")
    plt.grid()
    plt.ylim(0.5,8.5)
    plt.xlabel("Longitude")
    plt.title("MC12 Mixed Layer Depth")
    plt.colorbar(orientation="horizontal",ticks=np.arange(-8,11,4))
    plt.tight_layout()
    plt.show()
  elif plottype=="pcolor":
    # plot style 3: colour plot (this one was chosen)
    def smooth(data,n=36):
      x=convolve(data,np.ones((1,n))/n)
      return np.array([interp1d(hmix12.coord("longitude").points,x[i])(np.arange(87.5,161,5)) for i in range(8)])
    def smooth2(data,n=36,m1=False):
      x=convolve(data,np.ones((1,n))/n)
      if m1:
        return np.array([interp1d(hmix2.coord("longitude").points[1:-1],x[i])(np.arange(92.5,156,5)) for i in range(8)])
      return np.array([interp1d(hmix2.coord("longitude").points,x[i])(np.arange(92.5,156,5)) for i in range(8)])
    k=1
    lon = pcoloraxis(np.arange(87.5,161,5))
    lon2 = pcoloraxis(np.arange(92.5,156,5))
    ax=plt.subplot(241)
    plt.ylabel("MJO Phase")
    plt.pcolormesh(lon2,np.arange(0.5,9),smooth2(sst_r_mjo2)[:,::k],vmin=-0.225,vmax=0.225,cmap=ListedColormap(sns.color_palette("RdGy_r",11)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.xlabel("Longitude")
    plt.title("(a) MC2 Diurnal SST Range")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-0.2,0.21,0.1))
    ax=plt.subplot(242)
    a=plt.pcolormesh(lon2,np.arange(0.5,9),smooth2(sst_m_mjo2)[:,::k],vmin=-0.45,vmax=0.45,cmap=ListedColormap(sns.color_palette("bwr",11)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.title("(b) MC2 Foundational SST")
    plt.xlabel("Longitude")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-0.4,0.5,0.2))
    ax=plt.subplot(243)
    plt.pcolormesh(lon2,np.arange(0.5,9),smooth2(P_mjo2[:],m1=True)[:,::k],vmin=-9,vmax=9,cmap=ListedColormap(sns.color_palette("BrBG",11)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.xlabel("Longitude")
    plt.title("(c) MC2 Precip over Ocean")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-8,9,2))
    ax=plt.subplot(244)
    plt.pcolormesh(lon2,np.arange(0.5,9),smooth2(hmix_mjo2)[:,::k],vmin=-11,vmax=11,cmap=ListedColormap(sns.color_palette("PiYG",13)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.ylim(0.5,8.5)
    plt.xlabel("Longitude")
    plt.title("(d) MC2 Mixed Layer Depth")
    plt.colorbar(orientation="horizontal",ticks=np.arange(-8,11,4))
    ax=plt.subplot(245)
    plt.ylabel("MJO Phase")
    plt.pcolormesh(lon,np.arange(0.5,9),smooth(sst_r_mjo12)[:,::k],vmin=-0.225,vmax=0.225,cmap=ListedColormap(sns.color_palette("RdGy_r",11)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.xlabel("Longitude")
    plt.title("(e) MC12 Diurnal SST Range")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-0.2,0.21,0.1))
    ax=plt.subplot(246)
    a=plt.pcolormesh(lon,np.arange(0.5,9),smooth(sst_m_mjo12)[:,::k],vmin=-0.45,vmax=0.45,cmap=ListedColormap(sns.color_palette("bwr",11)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.title("(f) MC12 Foundational SST")
    plt.xlabel("Longitude")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-0.4,0.5,0.2))
    ax=plt.subplot(247)
    plt.pcolormesh(lon,np.arange(0.5,9),smooth(P_mjo12[:])[:,::k],vmin=-9,vmax=9,cmap=ListedColormap(sns.color_palette("BrBG",11)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.xlabel("Longitude")
    plt.title("(g) MC12 Precip over Ocean")
    plt.ylim(0.5,8.5)
    plt.colorbar(orientation="horizontal",ticks=np.arange(-8,9,2))
    ax=plt.subplot(248)
    plt.pcolormesh(lon,np.arange(0.5,9),smooth(hmix_mjo12)[:,::k],vmin=-11,vmax=11,cmap=ListedColormap(sns.color_palette("PiYG",13)))
    ax.set_yticks(np.arange(0.5,9),minor=True)
    ax.yaxis.grid(True,which="minor")
    ax.xaxis.grid(True)
    plt.ylim(0.5,8.5)
    plt.xlabel("Longitude")
    plt.title("(h) MC12 Mixed Layer Depth")
    plt.colorbar(orientation="horizontal",ticks=np.arange(-8,11,4))


    plt.tight_layout()
  plt.savefig(figname%plottype)
  plt.show()

#main("pcolor",True,"/work/scratch-pw2/emmah/eval/",figname="RMM_kpp_%s.png")
