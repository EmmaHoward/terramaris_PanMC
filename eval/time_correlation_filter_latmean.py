import iris
from iris.experimental.equalise_cubes import equalise_attributes
from iris.coord_categorisation import add_day_of_year
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
cy =iris.Constraint(latitude=lambda y: -15<=y<=15)
cx =iris.Constraint(longitude=lambda y: 85<=y<=155)
years = [2003,2014,2015,2016,2017]
path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/daily_latmean/"

bands = np.array([45,30,18,15,10,6,3])

def band_regress(data1,data2,bands):
  data1=detrend(data1,axis=1)
  data2=detrend(data2,axis=1)
  F1=np.fft.fft(data1,axis=1)                       
  F2=np.fft.fft(data2,axis=1)                       
  freq = np.fft.fftfreq(90)
  correl = []                      
  n=np.newaxis
#  fig,ax=plt.subplots(2,6) 
  for i,band in enumerate(bands):
    #mask = np.exp(-(np.abs(freq) - 1/band)**2/(0.025)**2)
    mask = np.abs(np.abs(freq) - 1/band)<0.005
#    print(mask)
    data3 = np.fft.ifft(F1*mask[np.newaxis,:,np.newaxis],axis=1).real
    data4 = np.fft.ifft(F2*mask[np.newaxis,:,np.newaxis],axis=1).real
#    ax[0,i].pcolormesh(data3[0])
#    ax[1,i].pcolormesh(data4[0])
    correl.append((((data3-data3.mean(axis=1)[:,n])*(data4-data4.mean(axis=1)[:,n])).mean(axis=1)[:,n]/data3.std(axis=1,ddof=0)[:,n]/data4.std(axis=1,ddof=0)[:,n]).mean())
    #correl.append( ((data3-data3.mean())*(data4-data4.mean())).mean()/data3.std(ddof=0)/data4.std(ddof=0))
  return np.array(correl)

def load_u(year,cy):
  ct = iris.Constraint(t = lambda tt: (tt.point.year,tt.point.month) in [(year,12),(year+1,1),(year+1,2)] and not (tt.point.month,tt.point.day) == (2,29))
  data1 = iris.load(["/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U200/U.%04d.global_domain.nc"%year,
                    "/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U200/U.%04d.global_domain.nc"%(year+1)],ct&cy&cx)
  data2 = iris.load(["/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U850/U.%04d.global_domain.nc"%year,
                    "/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U850/U.%04d.global_domain.nc"%(year+1)],ct&cy&cx)
  data1[1].coord("t").convert_units(data1[0].coord("t").units)
  data2[1].coord("t").convert_units(data2[0].coord("t").units)
  equalise_attributes(data1+data2)
  data1[0].coord("t").attributes={}
  data1[1].coord("t").attributes={}
  data2[0].coord("t").attributes={}
  data2[1].coord("t").attributes={}
  data1=data1.concatenate_cube()
  data2=data2.concatenate_cube()
  P200 = iris.coords.AuxCoord(200,units="hPa",standard_name="air_pressure")
  P850 = iris.coords.AuxCoord(850,units="hPa",standard_name="air_pressure")
  data1.add_aux_coord(P200)
  data2.add_aux_coord(P850)
  data = iris.cube.CubeList([data1,data2]).merge_cube()
  add_day_of_year(data,"t","doyr")
  print(data.coord("air_pressure"))
  data = data.aggregated_by("doyr",iris.analysis.MEAN)
  data = data.collapsed("latitude",iris.analysis.MEAN)
  return data



def main(var,i=None):
  MC12,obs = [],[]
  for year in years:
    cube = iris.load(path+"%04d_%s_15N15S.nc"%(year,var))[0]
    if var=="rainfall":
      cube = cube.rolling_window("longitude",iris.analysis.MEAN,7)
      ct = iris.Constraint(time = lambda t: (t.point.year,t.point.month) in [(year,12),(year+1,1),(year+1,2)] and not (t.point.month,t.point.day) == (2,29))
      gpm = iris.load(["/gws/nopw/j04/klingaman/datasets/GPM_IMERG/daily/3B-DAY.MS.MRG.3IMERG.%04d.V06.1x1.nc"%year,"/gws/nopw/j04/klingaman/datasets/GPM_IMERG/daily/3B-DAY.MS.MRG.3IMERG.%04d.V06.1x1.nc"%(year+1)],cx&cy&ct).concatenate_cube()
      obs.append(gpm.collapsed("latitude",iris.analysis.MEAN).data)
      lonobs = gpm.coord("longitude").points
      cube = cube.interpolate([("longitude",lonobs)],iris.analysis.Linear())
      MC12.append(cube.data)
    elif var=="olr":
        print("implement")
        exit()
    elif var=="eastward_wind":
        data = load_u(year,cy)
        lonobs = data.coord("longitude").points
        obs.append(data.data)
        cube = cube.rolling_window("longitude",iris.analysis.MEAN,5)
        cube = cube.interpolate([("longitude",lonobs)],iris.analysis.Linear())
        MC12.append(cube.data)
  if not i is None:
    MC12=np.array(MC12)[:,:90,i]
    obs=np.array(obs)[:,i]
  correl = band_regress(MC12,obs,bands)
  return correl


correl_u200 = main("eastward_wind",0)
correl_u850 = main("eastward_wind",1)
correl_P = main("rainfall")



def plot(correl):
    fig=plt.figure(figsize=(5,4))
    ax=plt.subplot(111)
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("Band-pass limited correlation")
    ax.set_ylim(-0.1,1)
    ax.set_xticks(1/bands[:])
    ax.grid()
    ax.set_xticklabels(["%d"%i for i in bands[:]],fontsize="x-small")
    shortname = ["U200","U850","precip"]
    for j in range(3):
      plt.plot(1/(bands[:]),correl[j],"o-",label=shortname[j])
    plt.legend()


plot([correl_u200,correl_u850,correl_P])
