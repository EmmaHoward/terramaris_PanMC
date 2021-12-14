import iris
import cftime
from iris.experimental.equalise_cubes import equalise_attributes
from panMC import panMC
from iris.coord_categorisation import add_day_of_year
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

cy =iris.Constraint(latitude=lambda y: -15<=y<=15)
cx =iris.Constraint(longitude=lambda y: 85<=y<=155)
cz = iris.Constraint(pressure=lambda p: p in [200,850])
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
    data3 = np.fft.ifft(F1*mask[np.newaxis,:,np.newaxis,np.newaxis],axis=1).real
    data4 = np.fft.ifft(F2*mask[np.newaxis,:,np.newaxis,np.newaxis],axis=1).real
#    ax[0,i].pcolormesh(data3[0])
#    ax[1,i].pcolormesh(data4[0])
    #correl.append((((data3-data3.mean(axis=1)[:,n])*(data4-data4.mean(axis=1)[:,n])).mean(axis=1)[:,n]/data3.std(axis=1,ddof=0)[:,n]/data4.std(axis=1,ddof=0)[:,n]).mean())
    correl.append( ((data3-data3.mean())*(data4-data4.mean())).mean()/data3.std(ddof=0)/data4.std(ddof=0))
  return np.array(correl)

def load_u(year,cy,dates):
  ct = iris.Constraint(t = lambda tt: (tt.point.year,tt.point.month) in [(year,12),(year+1,1),(year+1,2)] and not (tt.point.month,tt.point.day) == (2,29))
  data1 = iris.load(["/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U200/U.%04d.global_domain.nc"%year,
                    "/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U200/U.%04d.global_domain.nc"%(year+1)],cy&cx&ct)
  data2 = iris.load(["/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U850/U.%04d.global_domain.nc"%year,
                    "/gws/nopw/j04/klingaman/datasets/ERA-INTERIM/U850/U.%04d.global_domain.nc"%(year+1)],cy&cx&ct)
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
  print(data.coord("air_pressure"))
  return data


def read_data(var):
  MC12,obs = [],[]
  for year in years:
    if var=="rainfall":
      data = panMC(year,"MC12","rainfall").load_iris(Constraints=cy,variables=["convective_rainfall_amount","stratiform_rainfall_amount"])
      data = data[0]+data[1]
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
        ct = iris.Constraint(time=lambda t: t.point.hour in [0,6,12,18])
        cube = panMC(year,"MC12","eastward_wind_pd").load_iris(Constraints=cy&cz&ct)[0]
        add_day_of_year(cube,"time","doyr")
        cube = cube.aggregated_by("doyr",iris.analysis.MEAN)
        data = load_u(year,cy,cube.coord("time").units.num2date(cube.coord("time").points))
        obs.append(data.data)
        data.coord("longitude").guess_bounds()
        data.coord("latitude").guess_bounds()
        cube.coord("longitude").guess_bounds()
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").coord_system = data.coord("longitude").coord_system
        cube.coord("latitude").coord_system = data.coord("latitude").coord_system
        cube = cube.regrid(data,iris.analysis.AreaWeighted())
        MC12.append(cube.data)
  return MC12, obs

#correl_P = main("rainfall")



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
    for j in range(len(correl)):
      plt.plot(1/(bands[:]),correl[j],"o-",label=shortname[j])
    plt.legend()



u12,uobs = read_data("eastward_wind")
correl_u200 = band_regress(u12[:,:,0],uobs[:,0])
correl_u850 = band_regress(u12[:,:,1],uobs[:,1])

plot([correl_u200,correl_u850])#,correl_P])
