import iris
from scipy.fftpack import dct
import cftime
from iris.experimental.equalise_cubes import equalise_attributes
from panMC import panMC
from iris.coord_categorisation import add_day_of_year
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import iris.cube
cy =iris.Constraint(latitude=lambda y: -15<=y<=15)
cx =iris.Constraint(longitude=lambda y: 92<=y<152)
cz = iris.Constraint(pressure=lambda p: p in [200,850])
years = [2003,2014,2015,2016,2017,2018]
path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/daily_latmean/"

bands = np.array([45,30,18,15,10,6,3])


bins_x = np.array([150,90,50,35,27,22,18,16,14,12,10.5,9.5,8.5,7.5,6.5,5.5,4.5,3.5])
bins_t = np.array([135,55,   35,  25,  20,  16,  14,  12, 10.5, 9.5, 8.5,7.5,6.5,5.5,4.5,3.5,2.5,1.5])
tticks = np.array([30,15,10,5,2])
xticks = np.array([30,20,15,10,5])

def adjust_doyr(cube):
  doyr = cube.coord("doyr").points
  if (cube.coord("year").points%4==0).all():
    doyr[doyr>180] -= 1
  doyr[doyr<180] += 365
  cube.coord("doyr").points = doyr



def band_regress(data1,data2):
  correl=[0,0,0]
  for i,(y0,y1) in enumerate([(5,15),(-5,5),(-15,-5)]): 
    cy2 = iris.Constraint(latitude=lambda y: y0<=y<=y1)
    x1 = data1.extract(cy2).collapsed("latitude",iris.analysis.MEAN) 
    x2 = data2.extract(cy2).collapsed("latitude",iris.analysis.MEAN) 
    x1.data=detrend(x1.data.filled(0),axis=1)
    x2.data=detrend(x2.data.filled(0),axis=1)
    F1 = np.fft.fft(x1.data,axis=1)
    F2 = np.fft.fft(x2.data,axis=1)
    correl[i] = ((F1*F2.conj())/(np.abs(F1)*np.abs(F2))).mean(axis=0)
  correl=np.array(correl)
  return correl
     
def read_data(var):
  MC12,obs = iris.cube.CubeList(),iris.cube.CubeList()
  for year in years:
    print(year)
    yc = iris.coords.DimCoord(year,long_name="year",units="years")
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
        cube = panMC(year,"MC12","eastward_wind_pd").load_iris(Constraints=cz&ct)[0]
        add_day_of_year(cube,"time","doyr")
        cube = cube.aggregated_by("doyr",iris.analysis.MEAN)
#
        data = iris.load("/gws/nopw/j04/terramaris/emmah/era5/uv_%04d%02d.nc"%(year,(year+1)%100),cx&cy).extract(var)[0][:90]
        data = data.rolling_window("latitude",iris.analysis.MEAN,2)
        data.add_aux_coord(yc)
        add_day_of_year(data,"time","doyr")
        adjust_doyr(data)
        iris.util.promote_aux_coord_to_dim_coord(data,"doyr")
        data.remove_coord("time")
        data.rename("%s_obs"%var)   
        obs.append(data)
#
        data.coord("longitude").guess_bounds()
        cube.coord("longitude").guess_bounds()
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").coord_system = data.coord("longitude").coord_system
        cube.coord("latitude").coord_system = data.coord("latitude").coord_system
        cube = cube.regrid(data,iris.analysis.AreaWeighted())
        cube.add_aux_coord(yc)
        adjust_doyr(cube)
        cube.remove_coord("time")
        iris.util.promote_aux_coord_to_dim_coord(cube,"doyr")
        cube.remove_coord("forecast_reference_time")
        cube.rename("%s_MC12"%var)
        MC12.append(cube)
        iris.save([cube,data],"/work/scratch-pw2/emmah/%04d_winds_regridded_MC12.nc"%year)
  iris.util.equalise_attributes(obs)
  MC12 = MC12.merge_cube()
  obs = obs.merge_cube()
  return MC12, obs

#correl_P = main("rainfall")



def plot(correl,lon):#,correlx,correlt):
  fig = plt.figure()
  freqs = np.fft.fftfreq(90)
  for i in range(3):
    ax=plt.subplot(1,3,i+1)
    ax.grid()
    a=ax.contourf(lon,freqs[:45],correl[i][:45].real,cmap="PiYG",vmin=-1,vmax=1)
    plt.colorbar(a,orientation="horizontal")
  plt.show()

def read_data2(var):
  MC12,obs = iris.cube.CubeList(),iris.cube.CubeList()
  for year in years:
    yc = iris.coords.DimCoord(year,long_name="year",units="years")
    data=iris.load("/work/scratch-pw2/emmah/%04d_winds_regridded_MC12.nc"%year,cx)
    MC12.append(data.extract("eastward_wind_MC12")[0])
    obs.append(data.extract("eastward_wind_obs")[0])
  iris.util.equalise_attributes(obs)
  MC12 = MC12.merge_cube()
  obs = obs.merge_cube()
  return MC12, obs

import pickle
def main():
  u12,uobs = read_data2("eastward_wind")
  correl_u850 = band_regress(u12[:,:90,1],uobs[:,:,0])
#  pickle.dump((correl_u850,correlx_u850,correlt_u850),open("correl_u_850.pickle","wb"))
#  correl_u850,correlx_u850,correlt_u850 = pickle.load("correl_u_850.pickle","rb")
  plot(correl_u850,uobs.coord("longitude").points)
#  correl_u200,correlx_u200,correlt_u200 = band_regress(u12[:,:90,0],uobs[:,:,1])
#  pickle.dump((correl_u200,correlx_u200,correlt_u200),open("correl_u_200.pickle","wb"))
#  plot(correl_u200,correlx_u200,correlt_u200)
  plt.show()
  import pdb;pdb.set_trace()

main()
