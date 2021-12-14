#!/home/users/emmah/miniconda3/envs/iris3/bin/python
#SBATCH -p short-serial
#SBATCH -o /home/users/emmah/log/corr-%a.o
#SBATCH -e /home/users/emmah/log/corr-%a.e 
#SBATCH --array=[1-90]
#SBATCH -t 02:00:00
from scipy.fftpack import dct
from matplotlib.cm import get_cmap
import iris
import numpy as np
import datetime as dt
from panMC import panMC
from iris.experimental.equalise_cubes import equalise_attributes
cp = iris.Constraint(air_pressure = lambda p: p in [850,500,200])
cx = iris.Constraint(longitude=lambda lon: 92.5<=lon<152.5)
cy = iris.Constraint(latitude=lambda lat: -15<=lat<=15)
from scipy.ndimage import convolve
import matplotlib.pyplot as plt
bands = np.array([30,15,10,7.5,6,5,3.75,3,2.5,2,1.5,1])#,0.5,0.25]
from iris.coord_categorisation import add_day_of_year
import os
from scipy.signal import detrend


def pad(var,data,buf):
  out = iris.cube.CubeList()
  for cube in data:
    cube.data=detrend(cube.data.filled(0),axis=0)
    lat = cube.coord("latitude").points
    y0 = lat.min() + buf
    y1 = lat.max() - buf
    lon = cube.coord("longitude").points
    x0 = lon[0] + buf
    x1 = lon[-1] - buf
    s1 = (np.exp(-( ((lon>x1)*(lon-x1) - (lon<x0)*(lon-x0))**2 )/buf**2))
    s2 = (np.exp(-( ((lat>y1)*(lat-y1) - (lat<y0)*(lat-y0))**2 )/buf**2))
    s12  = np.meshgrid(s1,s2)
    support = s12[0]*s12[1]
#    support = iris.coords.AuxCoord(s12[0]*s12[1],units=1)
    if var == "rainfall":
      out.append(iris.analysis.maths.multiply(cube,support,(1,2)))
    else: 
      anom = cube - cube.collapsed(["time"],iris.analysis.MEAN)
      out.append(iris.analysis.maths.multiply(anom,support,(2,3)))
  return out
  


def dct_correl(data,bands,dx=0.25,dy=0.25,pl = True,axes="xyt"):
  transformed = np.array([np.fft.fft2(cube.data.filled(0)) for cube in data])
  ny,nx = len(data[0].coord("latitude").points),len(data[0].coord("longitude").points)
  freqy = np.fft.fftfreq(ny,dy)
  freqx = np.fft.fftfreq(nx,dx)
  freqxx,freqyy=np.meshgrid(freqx,freqy)
  k = np.abs(np.hypot(freqxx,freqyy))
  out = []
  n=np.newaxis
  for i in range(len(bands[:])):
    k0 = 1/bands[i]
    mask = np.exp(-( k-k0)**2/(0.5/ny/dy)**2)
    filtered = np.fft.ifft2(transformed*mask)
    #grid point correl
    if axes=="t":
      out.append((((filtered - filtered.mean(axis=1)[:,np.newaxis])[:,np.newaxis]*(filtered - filtered.mean(axis=1)[:,np.newaxis])[np.newaxis]).mean(axis=2)
                       /filtered.std(axis=1,ddof=0)[:,np.newaxis]/filtered.std(axis=1,ddof=0)[np.newaxis]).mean(axis=(-1,-2)))
    # x y t correl
    elif axes=="xyt":
      out.append(((filtered - filtered.mean(axis=(1,-2,-1))[:,n,...,n,n])[:,n]*(filtered - filtered.mean(axis=(1,-2,-1))[:,n,...,n,n])[n]).mean(axis=(2,-2,-1))
                       /filtered.std(axis=(1,-2,-1),ddof=0)[:,n]/filtered.std(axis=(1,-2,-1),ddof=0)[n])
    c = list(get_cmap("viridis")(i/len(bands)))
  out = np.array(out)
  plt.show()
  return out

def to_daily(data):
  out = iris.cube.CubeList()
  for cube in data:
    add_day_of_year(cube,"time","doyr")
    out.append(cube.aggregated_by("doyr",iris.analysis.MEAN))
  return out

def regrid(data):
  data2,data12,era5=data 
  template = era5[:,:,1:].copy()
  template.coord("latitude").points = ((era5.coord("latitude").points[1:]+era5.coord("latitude").points[:-1]))/2
  era5 = era5.regrid(template,iris.analysis.Linear())
  era5.coord("longitude").coord_system=data12.coord("longitude").coord_system
  era5.coord("latitude").coord_system=data12.coord("latitude").coord_system
  era5.coord("longitude").guess_bounds()
  era5.coord("latitude").guess_bounds()
  era5=era5[:,::-1]
  data12.coord("longitude").guess_bounds()
  data12.coord("latitude").guess_bounds()
  data2.coord("longitude").guess_bounds()
  data2.coord("latitude").guess_bounds()
  data2 = data2.regrid(era5,iris.analysis.AreaWeighted())
  data12 = data12.regrid(era5,iris.analysis.AreaWeighted())
  import pdb;pdb.set_trace()
  return (data2,data12,era5)

def load(var,pressure):
  cp =iris.Constraint(pressure=lambda p: p in pressure)
  cp2 =iris.Constraint(pressure_level=lambda p: p in pressure)
  dates = [dt.datetime(2015,12,1)+dt.timedelta(i) for i in range(31)]
  data12 = panMC(2015,"MC12",var+"_pd").load_iris(dates,Constraints=cp)[0]
  data2 = panMC(2015,"MC2-tmp",var+"_pd").load_iris(dates,Constraints=cp)[0]
  ct = iris.Constraint(time = lambda t: dt.datetime(t.point.year,t.point.month,t.point.day) in dates)
  era5 = iris.load(["/gws/nopw/j04/terramaris/emmah/era5/uvq_%04d%02d.nc"%(y,m) for (y,m) in [(2015,12),(2016,1)]],cx&cy&ct&cp2).extract(var)[0]
  data2.coord("time").points = data2.coord("time").points - 0.5
  data12.coord("time").points = data12.coord("time").points - 0.5
  return (data2,data12,era5)

def main(var,pressure,axes):
  if os.path.exists("/work/scratch-nopw/emmah/%s_201512_0p25.nc"%var):
    data=iris.load("/work/scratch-nopw/emmah/%s_201512_0p25.nc"%var)
    data = data.extract([var+"_MC2",var+"_MC12",var+"_era5"])
  else:
    data = load(var,pressure)
    print("loaded")
    data = to_daily(data)
    print("daily means computed")
    data = regrid(data)
    print("regridded")
    print(data)
    data[0].rename(var+"_MC2")
    data[1].rename(var+"_MC12")
    data[2].rename(var+"_era5")
    iris.save(data,"/work/scratch-nopw/emmah/%s_201512_0p25.nc"%var,zlib=True)
  data=pad(var,data,5)
  correl = dct_correl(data,bands,axes=axes)
  print("computed")
  plot(correl,pressure,var,axes)

def plot(correl,pressure,var,axes):
  fig=plt.figure(figsize=(8,4))
  for i in range(2):
    ax=plt.subplot(1,2,i+1)
    ax2 =ax.twiny()
    ax.set_xlabel("Zonal Wavenumber")
    ax2.set_xlabel("Wavelength (degrees)")
    if i==0:
      ax.set_ylabel("Band-pass limited correlation")
    ax.set_ylim(-0.1,1)
    ax.set_xlim(0,62)
    ax2.set_xlim(0,62)
    ax2.set_xticks(60/bands[1::2])
    ax.set_xticks(60/bands[1::2])
    ax.grid()
    ax2.set_xticklabels(bands[1::2],fontsize="x-small")
    shortname = {"northward_wind":"V","eastward_wind":"U","specific_humidity":"Q"}[var]
    for j in range(3):
      i1,i2 = [[0,1],[0,2],[1,2]][j]
      name = ["MC2 and MC12","MC2 and ERA5","MC12 and ERA5"][j]
      plt.plot(60/(bands[:]),correl[:,i1,i2,i],"o-",label=name)
    plt.legend(title="%s%d"%(shortname,pressure[i]))
  fig.savefig("figs/correl_%s_%s.png"%(shortname,axes))

#main("eastward_wind",[300,850],"xyt")
#main("eastward_wind",[300,850],"t")
#main("northward_wind",[300,850],"xyt")
main("northward_wind",[300,850],"t")


plt.show()
