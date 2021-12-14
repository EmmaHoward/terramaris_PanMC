from iris.util import equalise_attributes
import iris
from panMC import panMC
import os
from iris.coord_categorisation import add_day_of_year
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
from bias_plots import choose_bounds
cy =iris.Constraint(latitude=lambda y: -15<=y<=15)
cx =iris.Constraint(longitude=lambda y: 85<=y<=160)

def adjust_doyr(cube):
  doyr = cube.coord("doyr").points
  if (cube.coord("year").points%4==0).all():
    doyr[doyr>180] -= 1
  doyr[doyr<180] += 365
  cube.coord("doyr").points = doyr
obsname = {"rainfall":"GPM-IMERG","eastward_wind":"ERA-Interim","olr":"ERA-5"}
MC2_path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/"
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/"

def latmeanMC12(year,var,stream,Constraints,name):
  if os.path.exists(MC12_path+"daily_latmean/"+"%04d_%s_15N15S.nc"%(year,name)):
    data = iris.load(MC12_path+"daily_latmean/"+"%04d_%s_15N15S.nc"%(year,name))[0]
    return data
  data = panMC(year,"MC12",stream).load_iris(Constraints=Constraints,variables=var)
  print(data)
  data=data[0].collapsed("latitude",iris.analysis.MEAN)
  add_day_of_year(data,"time","doyr")
  data = data.aggregated_by("doyr",iris.analysis.MEAN)
  iris.save(data,path+"%04d_%s_15N15S.nc"%(year,name),zlib=True)
  return data

  c_obs,base = lagged_regress(obs,lons,ndays,base=None,sign=sign,lowpass=lowpass,cutoff=cutoff)

def latmeanMC2(year,var,stream,Constraints,name):
  if 1:#name=="rainfall":
    data=iris.load([MC2_path+"/precip/%04d%02d_hourly_coarsened_rainfall.nc"%(year+int(month<6),month) for month in [12,1,2] ])
    equalise_attributes(data)
    data = data.concatenate_cube()*4
  elif name=="eastward_wind":
    print("still running")
    exit()
  else:
    print("not yet implemented")
    exit()
  data = data.collapsed("latitude",iris.analysis.MEAN).aggregated_by("doyr",iris.analysis.MEAN)
  return data


def lagged_regress(data,lons,ndays,sign=1,basedata=None,lowpass=True,cutoff=20):
  data2=detrend(data.data,axis=1)
  if data2.ndim==2:
    data2=np.array([data2])
  print(data2.shape)
  F=np.fft.fft2(data2,axes=(1,2))                       
  freqt = np.fft.fftfreq(90)                       
  freqx = np.fft.fftfreq(data2.shape[2])                       
  F2=F.copy()
  freqxx,freqtt = np.meshgrid(freqx,freqt)
  if lowpass:
    F2[:,np.abs(freqt)>1/cutoff]=0                       
  else:
    F2[:,np.abs(freqt)<1/cutoff]=0                       
  if sign != 0:
    F2 = F2*((freqtt<0)*((freqxx*sign)>0)+(freqtt>0)*((freqxx*sign)<0))
  data3=np.fft.ifft2(F2,axes=(1,2)).real                       
  correl = np.zeros((len(lons),2*ndays+1,data2.shape[2]))
  data3 = data3.transpose([1,0,2])
  save_basedata = []
  for j,lon in enumerate(lons):
    i_base = np.argmin(np.abs(data.coord("longitude").points-lon))
    if basedata is None:
      base = data3[:,:,i_base]
    else:
      base = basedata[j]
    save_basedata.append(base)
    correl[j,ndays] = (((data3-data3.mean(axis=0))*(base-base.mean(axis=0))[:,:,np.newaxis]).mean(axis=0)/data3.std(axis=0,ddof=0)/base.std(axis=0,ddof=0)[:,np.newaxis]).mean(axis=0)
    for i in range(1,ndays):
      correl[j,ndays-i] =  (((data3[:-i]-data3[:-i].mean(axis=0))*(base[i:]-base[i:].mean(axis=0))[:,:,np.newaxis]).mean(axis=0)/data3[:-i].std(axis=0,ddof=0)/base[i:].std(axis=0,ddof=0)[:,np.newaxis]).mean(axis=0)
      correl[j,ndays+i] =  (((data3[i:]-data3[i:].mean(axis=0))*(base[:-i]-base[:-i].mean(axis=0))[:,:,np.newaxis]).mean(axis=0)/data3[i:].std(axis=0,ddof=0)/base[:-i].std(axis=0,ddof=0)[:,np.newaxis]).mean(axis=0)
  save_basedata = np.array(save_basedata)
  return correl,save_basedata


def lagged_regress1D(data,i_base,n,freqlim,base=None):
  data2=detrend(data,axis=1)
  F=np.fft.fft(data2,axis=1)                       
  freq = np.fft.fftfreq(90)                       
  F2=F.copy()                       
  F2[:,np.abs(freq)>1/freqlim]=0                       
  data3=np.fft.ifft(F2,axis=1).real                       
  data3 = data3.transpose([1,0,2])
  if base is None:
    base = data3[:,:,i_base]
  correl = np.zeros((2*n+1,data.shape[2]))
  correl[n] = (((data3-data3.mean(axis=0))*(base-base.mean(axis=0))[:,:,np.newaxis]).mean(axis=0)/data3.std(axis=0,ddof=0)/base.std(axis=0,ddof=0)[:,np.newaxis]).mean(axis=0)
  for i in range(1,n):
    correl[n-i] =  (((data3[:-i]-data3[:-i].mean(axis=0))*(base[i:]-base[i:].mean(axis=0))[:,:,np.newaxis]).mean(axis=0)/data3[:-i].std(axis=0,ddof=0)/base[i:].std(axis=0,ddof=0)[:,np.newaxis]).mean(axis=0)
    correl[n+i] =  (((data3[i:]-data3[i:].mean(axis=0))*(base[:-i]-base[:-i].mean(axis=0))[:,:,np.newaxis]).mean(axis=0)/data3[i:].std(axis=0,ddof=0)/base[:-i].std(axis=0,ddof=0)[:,np.newaxis]).mean(axis=0)
  return correl,base

def load(var,years):
  MC2,MC12,obs = iris.cube.CubeList(),iris.cube.CubeList(),iris.cube.CubeList()
  for year in years:
    yc = iris.coords.DimCoord(year,long_name="year",units="years")
    if 0:#os.path.exists(path+"%04d_%s_15N15S.nc"%(year,var)):
    #cube = iris.load(path+"%04d_olr_15N15S.nc"%year)[0]
      cube = iris.load(path+"%04d_%s_15N15S.nc"%(year,var))[0]
    else:
      if var=="rainfall":
        cube12 = latmeanMC12(year,["convective_rainfall_amount","stratiform_rainfall_amount"],"rainfall",cy,var)
        cube2 = latmeanMC2(year,["stratiform_rainfall_amount"],"rainfall",cy,var)
      elif var == "eastward_wind":
        cz = iris.Constraint(pressure = lambda p: p in [850,200])
        cube12 = latmeanMC12(year,["eastward_wind"],"eastward_wind_pd",cy&cz,var)
        cube2 = latmeanMC2(year,["eastward_wind"],"eastward_wind_pd",cy&cz,var)
    cube12.add_aux_coord(yc)
    cube2.add_aux_coord(yc)
    adjust_doyr(cube2)
    adjust_doyr(cube12)
    iris.util.promote_aux_coord_to_dim_coord(cube12,"doyr")
    iris.util.promote_aux_coord_to_dim_coord(cube2,"doyr")
    cube2.remove_coord("time")
    cube12.remove_coord("time")
    if var=="rainfall":
      cube12 = cube12.rolling_window("longitude",iris.analysis.MEAN,7)
      cube2 = cube2.rolling_window("longitude",iris.analysis.MEAN,7)
      ct = iris.Constraint(time = lambda t: (t.point.year,t.point.month) in [(year,12),(year+1,1),(year+1,2)] and not (t.point.month,t.point.day) == (2,29))
      gpm = iris.load(["/gws/nopw/j04/klingaman/datasets/GPM_IMERG/daily/3B-DAY.MS.MRG.3IMERG.%04d.V06.1x1.nc"%year,"/gws/nopw/j04/klingaman/datasets/GPM_IMERG/daily/3B-DAY.MS.MRG.3IMERG.%04d.V06.1x1.nc"%(year+1)],cx&cy&ct).concatenate_cube()
      add_day_of_year(gpm,"time","doyr")
      gpm.add_aux_coord(yc)
      adjust_doyr(gpm)
      iris.util.promote_aux_coord_to_dim_coord(gpm,"doyr")
      gpm.remove_coord("time")
      obs.append(gpm.collapsed("latitude",iris.analysis.MEAN))
    elif var=="eastward_wind":
        data = iris.load("/gws/nopw/j04/terramaris/emmah/era5/uv_%04d%02d.nc"%(year,(year+1)%100),cx&cy).extract(var)[0][:90]
        data.add_aux_coord(yc)
        add_day_of_year(data,"time","doyr")
        adjust_doyr(data)
        iris.util.promote_aux_coord_to_dim_coord(data,"doyr")
        data.remove_coord("time")
        data = data.collapsed("latitude",iris.analysis.MEAN)
        obs.append(data)
    MC12.append(cube12)
    MC2.append(cube2)
  for cube in MC12:
    cube.remove_coord("forecast_reference_time")
  iris.util.equalise_attributes(obs)
  MC12=MC12.merge_cube()
  MC2=MC2.merge_cube()
  obs = obs.merge_cube()
  return MC12,MC2,obs



def pad(var,cube,buf):
  lon = cube.coord("longitude").points
  x0 = lon[0] + buf
  x1 = lon[-1] - buf
  support = iris.coords.AuxCoord(np.exp(-( ((lon>x1)*(lon-x1) - (lon<x0)*(lon-x0))**2 )/buf**2),units=1)
  if var == "rainfall":
    return  iris.analysis.maths.multiply(cube,support,1)
  else: 
    anom = cube - cube.collapsed(["year","doyr"],iris.analysis.MEAN)
    return  iris.analysis.maths.multiply(anom,support,1)
  

def plot(data,data2,ndays,longitudes,longitudes2,basepoints,names,titles,figname=None):
  fig,axs = plt.subplots(len(data[0]),len(data),figsize=(9,10))
  fig.subplots_adjust(top=0.95,bottom=0.15,left=0.09,right=0.969,hspace=0.189,wspace=0.305)
  for j,correl in enumerate(data):
    correl2 = data2[j]
    ticks = choose_bounds( np.percentile(correl[0],99),np.percentile(correl[0],1),12)
    for i in range(len(correl)):
      ax = axs[i,j]
#      a=ax.contour(longitudes2[i],range(-1*ndays[j],ndays[j]+1),correl2[i][0],np.arange(-0.875,1,0.25),colors=["k"],linewidths=[1],extend="both")
      a=ax.contourf(longitudes[i],range(-1*ndays[j],ndays[j]+1),correl[i][0],ticks,cmap="PiYG",extend="both")
      ax.plot([basepoints[j]],[0],".c")
      ax.grid(True)
      if j==0:
        ax.set_ylabel("%s \n Lead (days)      Lag (days)"%names[i],fontsize="small")
      if i==0:
        ax.set_title(titles[j])
      if i==len(correl)-1:
        ax.set_xlabel("Longitude")
    cax=fig.add_axes([ax.get_position().x0,0.04,ax.get_position().width,0.05])
    cb = fig.colorbar(a,cax=cax,orientation="horizontal",ticks = 0.5*(ticks[1:]+ticks[:-1])[::len(ticks)//5])
  if not figname is None:
    fig.savefig(figname)


def main(years,figname):
  """
  var="eastward_wind"
  MC12,obs = load(var,years)
  print(MC12)
  print(obs)
  MC12 = pad(var,MC12[:,:-1,1],5)
  obs = pad(var,obs[:,:,0],5)
  w_lonObs,w_lon12 = [obs.coord("longitude").points,MC12.coord("longitude").points]

  wind_MJO_obs,wind_MJO_base = lagged_regress(obs,[100],25,basedata=None,         sign=0,lowpass=True, cutoff=20)
  wind_MJO_MC12_bpm          = lagged_regress(MC12,[100],25,basedata=None,         sign=0,lowpass=True, cutoff=20)[0] # base point model
  wind_MJO_MC12_bpo          = lagged_regress(MC12,[100],25,basedata=wind_MJO_base,sign=0,lowpass=True, cutoff=20)[0] # base point obs
  wind_MJO = [wind_MJO_obs,wind_MJO_MC12_bpm,wind_MJO_MC12_bpo]

  wind_Kel_obs,wind_Kel_base = lagged_regress(obs,[100],15,basedata=None,         sign=1,lowpass=False,cutoff=20)
  wind_Kel_MC12_bpm =          lagged_regress(MC12,[100],15,basedata=None,         sign=1,lowpass=False,cutoff=20)[0] # base point model
  wind_Kel_MC12_bpo =          lagged_regress(MC12,[100],15,basedata=wind_Kel_base,sign=1,lowpass=False,cutoff=20)[0] # base point obs
  wind_Kel =  [wind_Kel_obs,wind_Kel_MC12_bpm,wind_Kel_MC12_bpo]


  wind_Ros_obs,wind_Ros_base = lagged_regress(obs,[140],15,basedata=None,         sign=-1,lowpass=False,cutoff=20)
  wind_Ros_MC12_bpm =          lagged_regress(MC12,[140],15,basedata=None,         sign=-1,lowpass=False,cutoff=20)[0] # base point model
  wind_Ros_MC12_bpo =          lagged_regress(MC12,[140],15,basedata=wind_Ros_base,sign=-1,lowpass=False,cutoff=20)[0] # base point obs
  wind_Ros =  [wind_Ros_obs,wind_Ros_MC12_bpm,wind_Ros_MC12_bpo]
  """
  var="rainfall"
  MC12,MC2,obs = load(var,years)
  MC12 = pad(var,MC12[:,:],5)
  MC2 = pad(var,MC2[:,:],5)
  obs = pad(var,obs[:],5)
  r_lonObs,r_lon12,r_lon2 = [obs.coord("longitude").points,MC12.coord("longitude").points,MC2.coord("longitude").points]

  rain_MJO_obs,rain_MJO_base = lagged_regress(obs,[100],25,basedata=None,         sign=0,lowpass=True, cutoff=20)
  rain_MJO_MC12_bpm          = lagged_regress(MC12,[100],25,basedata=None,         sign=0,lowpass=True, cutoff=20)[0] # base point model
  rain_MJO_MC12_bpo          = lagged_regress(MC12,[100],25,basedata=rain_MJO_base,sign=0,lowpass=True, cutoff=20)[0] # base point obs
  rain_MJO_MC2_bpm          = lagged_regress(MC2,[100],25,basedata=None,         sign=0,lowpass=True, cutoff=20)[0] # base point model
  rain_MJO_MC2_bpo          = lagged_regress(MC2,[100],25,basedata=rain_MJO_base,sign=0,lowpass=True, cutoff=20)[0] # base point obs
  rain_MJO = [rain_MJO_obs,rain_MJO_MC2_bpm,rain_MJO_MC2_bpo,rain_MJO_MC12_bpm,rain_MJO_MC12_bpo]

  rain_Kel_obs,rain_Kel_base = lagged_regress(obs,[100],15,basedata=None,         sign=1,lowpass=False,cutoff=20)
  rain_Kel_MC12_bpm =          lagged_regress(MC12,[100],15,basedata=None,         sign=1,lowpass=False,cutoff=20)[0] # base point model
  rain_Kel_MC12_bpo =          lagged_regress(MC12,[100],15,basedata=rain_Kel_base,sign=1,lowpass=False,cutoff=20)[0] # base point obs
  rain_Kel_MC2_bpm =          lagged_regress(MC2,[100],15,basedata=None,         sign=1,lowpass=False,cutoff=20)[0] # base point model
  rain_Kel_MC2_bpo =          lagged_regress(MC2,[100],15,basedata=rain_Kel_base,sign=1,lowpass=False,cutoff=20)[0] # base point obs
  rain_Kel = [rain_Kel_obs,rain_Kel_MC2_bpm,rain_Kel_MC2_bpo,rain_Kel_MC12_bpm,rain_Kel_MC12_bpo]

  rain_Ros_obs,rain_Ros_base = lagged_regress(obs,[140],15,basedata=None,         sign=-1,lowpass=False,cutoff=20)
  rain_Ros_MC12_bpm =          lagged_regress(MC12,[140],15,basedata=None,         sign=-1,lowpass=False,cutoff=20)[0] # base point model
  rain_Ros_MC12_bpo =          lagged_regress(MC12,[140],15,basedata=rain_Ros_base,sign=-1,lowpass=False,cutoff=20)[0] # base point obs
  rain_Ros_MC2_bpm =          lagged_regress(MC2,[140],15,basedata=None,         sign=-1,lowpass=False,cutoff=20)[0] # base point model
  rain_Ros_MC2_bpo =          lagged_regress(MC2,[140],15,basedata=rain_Ros_base,sign=-1,lowpass=False,cutoff=20)[0] # base point obs
  rain_Ros = [rain_Ros_obs,rain_Ros_MC2_bpm,rain_Ros_MC2_bpo,rain_Ros_MC12_bpm,rain_Ros_MC12_bpo]

  wind_MJO,wind_Kel,wind_Ros = rain_MJO,rain_Kel,rain_Ros
  w_lonObs,w_lon12,w_lon2 = r_lonObs,r_lon12,r_lon2
  
  plot([rain_MJO,rain_Kel,rain_Ros],[wind_MJO,wind_Kel,wind_Ros],[25,15,15],[r_lonObs,r_lon2,r_lon2,r_lon12,r_lon12],[w_lonObs,w_lon2,w_lon2,w_lon12,w_lon12],[100,100,140],["GPM & ERA-5","MC2","MC2 vs Obs","MC12","MC12 vs Obs"],["MJO \n (>20 days)","Kelvin \n (<20 days, eastward)","Rossby \n (<20 days, westward)"],figname)

if __name__=="__main__":
  main([2016])
