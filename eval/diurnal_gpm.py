#
# plot rainfall diurnal cycles (Figure 5)
#

import iris
import iris.plot as iplt
from panMC import panMC
import cmocean
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import subprocess
from iris.coord_categorisation import add_hour,add_month
from iris.util import equalise_attributes
from iris.coord_categorisation import add_categorised_coord
from subprocess import check_call
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import sst_bias
import precip_bias
from matplotlib.cm import get_cmap
from scipy.interpolate import CubicSpline


# domain limits
cx = iris.Constraint(longitude = lambda x: 84 < x < 161)
cy = iris.Constraint(latitude  =  lambda y: -21 < y < 21)


MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/precip/"
MC2_path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"
ref_path = "/gws/nopw/j04/klingaman/datasets/GPM_IMERG/monthly/"

from iris.coord_categorisation import add_hour,add_month

def add_time_of_day(cube, coord, name='hour'):
    def _time_of_day(coord, value):
        pt = coord.units.num2date(value)
        return pt.hour+pt.minute/60.0+pt.second/3600
    add_categorised_coord(cube, name, coord,
                          _time_of_day)
	
def diurnal_amp_peak(cube):
  # calculate time of diurnal max in local standard time, diurnal range and overall mean precip
  cube = cube.aggregated_by("hour",iris.analysis.MEAN)
  lon = cube.coord("longitude").points
  utc = cube.coord("hour").points
#    tmp = np.meshgrid(lon,utc)
#    lst = tmp[1]+tmp[0]/180*12
  lst = (utc[cube.data.argmax(axis=0)]+lon/180*12)%24
  lst =  cube[0].copy(data=lst)
  lst.units="hours"
  mean = cube.collapsed("time",iris.analysis.MEAN)
  amp = (cube.collapsed("time",iris.analysis.MAX) - cube.collapsed("time",iris.analysis.MIN))/mean
  return amp,lst,mean

def load_precip(MC12_years,MC2_years):
  # load rainfall data (use precip_bias loader)
  P12 = iris.cube.CubeList()
  P2 = iris.cube.CubeList()
  for year in MC12_years:
    print(year)
    P12.append(precip_bias.load(year,"MC12",MC12_path))
    P12[-1].coord("time").convert_units("days since 2003-01-01") 
  P12=P12.concatenate_cube()
  for year in MC2_years:
    print(year)
    #P2.append(iris.load(MC2_path%(year,(year+1)%100))[0])
    P2.append(precip_bias.load(year,"MC2",MC2_path))
    P2[-1].coord("time").convert_units("days since 2003-01-01") 
    #P2[-1].cell_methods = ()
  iris.util.equalise_attributes(P2)
  P2=P2.concatenate_cube()
  data=P2.data
  #data[24:]*=15
  P2.data=data
  return P2,P12


def plot(years12,years2,figname):
  #plot diurnal cycle figures
  # load data
  P2,P12 = load_precip(years12,years2)
  P2 = P2.regrid(P12,iris.analysis.AreaWeighted())
  GPM = iris.load_cube("/work/scratch-pw2/emmah/gpm_diurnal.nc")
  GPM.data = GPM.data[:]
  GPM.convert_units(P2.units)
  # compute diurnal cycle properties
  amp12,peak12,mean12 = diurnal_amp_peak(P12)
  amp2,peak2,mean2 = diurnal_amp_peak(P2)
  ampref,peakref,meanref = diurnal_amp_peak(GPM)
  # mask data with weak diurnal cycle or low rain rates
  peak_2  = peak2.copy()
  peak_2.data =np.ma.masked_array(peak_2.data,(amp2.data<1)+(mean2.data<5))
  peak_12  = peak12.copy()
  peak_12.data =np.ma.masked_array(peak_12.data,(amp12.data<1)+(mean12.data<5))
  peak_ref = peakref.copy()
  peak_ref.data = np.ma.masked_array(peakref.data,(ampref.data<1)+(meanref.data<5))
  # build colormap, adjust to make more cyclic
  rgb = get_cmap('turbo')(list(np.arange(1/24,1,1/12)))
  rgb_new = rgb.copy()
  n,skip = 12,2
  for i in range(3):
    rgb_new[:,i] = CubicSpline(list(range(0,n//2-skip))+list(range(n//2+skip,n+1)),rgb[list(range(n//2,n-skip))+list(range(skip,n//2+1))][:,i])(list(range(n//2,n))+list(range(0,n//2)))
  cdict = {colour: np.array([np.linspace(0,1,13),[rgb_new[0,i]*0]+list(rgb_new[:,i]),list(rgb_new[:,i])+[rgb_new[-1,i,]]]).T for i,colour in enumerate(['red','green','blue'])}
  cmap1=LinearSegmentedColormap("",cdict)
  fig=plt.figure(figsize=(10,3))
  # build figure
  ax,p = [],[]
  for i,(cube,cmap,vmin,vmax,title,extend) in enumerate([(peak_2    ,cmap1,  0,24  ,"(a) MC2 Hour of Max Precipitation",None),
                                     (peak_12    ,cmap1,  0,24,  "(b) MC12 Hour of Max Precipitation",None),
                                     (peak_ref   ,cmap1,  0,24,  "(c) GPM Hour of Max Precipitation",None)]):
    ax.append(plt.subplot(1,3,1+i,projection=ccrs.PlateCarree()))
    ax[-1].coastlines()
    plt.xlim(90,155)
    plt.ylim(-15,15)
    p.append(iplt.pcolormesh(cube,vmin=vmin,vmax=vmax,cmap=cmap))
    plt.title(title)
  fig.colorbar(p[0],ax=ax,ticks=np.arange(0,25,4),orientation="horizontal")
  fig.subplots_adjust(top=0.95,bottom=0.35,left=0.015,right=0.985,hspace=0.2,wspace=0.048)
  fig.savefig(figname,dpi=300)
  plt.show()

years= [2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]
if __name__ == "__main__":
  plot(years,years,"/home/users/emmah/eval_figs_Feb23/diurnal_cycle_new.png")


