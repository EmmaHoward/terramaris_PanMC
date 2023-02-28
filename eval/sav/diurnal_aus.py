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
from matplotlib.colors import ListedColormap
import sst_bias
import precip_bias
 
cx = iris.Constraint(longitude = lambda x: 84 < x < 161)
cy = iris.Constraint(latitude  =  lambda y: -21 < y < 21)


#MC12_years = [2003,2014,2015,2016,2017]
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/precip/"
MC2_path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"

#MC2_years = [2015,2016]
#MC2_path = "/gws/nopw/j04/terramaris/emmah/Q1Q2_analysis/%04d%02d_MC2_coarsened_diurnal_precip_MC12.nc" 

ref_path = "/gws/nopw/j04/klingaman/datasets/GPM_IMERG/monthly/"


def add_time_of_day(cube, coord, name='hour'):
    def _time_of_day(coord, value):
        pt = coord.units.num2date(value)
        return pt.hour+pt.minute/60.0+pt.second/3600
    add_categorised_coord(cube, name, coord,
                          _time_of_day)
	
def diurnal_amp_peak(cube):
#  add_time_of_day(cube,"time","minhour")
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


def load_sst(years12,years2):
#  sst2 = iris.cube.CubeList()
#  sst2.append(sst_bias.load(2015,"MC2-tmp","/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/postprocessed_outputs/"))
#  sst2.append(sst_bias.load(2016,"MC2","/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/"))
#  iris.util.equalise_attributes(sst2)
#  for cube in sst2:
#    cube.coord("time").convert_units("days since 2003-01-01")
#  sst2=sst2.concatenate_cube()

  sst12,sst2,sstref = sst_bias.load_all(years12,years2)
#  sst12 = sst12.extract(iris.Constraint(time=lambda t: t.point.month==12))
#  print(sst2)
#  print(sst12)
  sst2 = sst2[:,0].aggregated_by("hour",iris.analysis.MEAN)
  sst12 = sst12[:,0].aggregated_by("hour",iris.analysis.MEAN)
  return sst2,sst12

def plot(years12,years2,figname):
  P2,P12 = load_precip(years12,years2)
  P2 = P2.regrid(P12,iris.analysis.AreaWeighted())
  sst2,sst12 = load_sst(years12,years2)
  sst2_range = sst2.collapsed("hour",iris.analysis.MAX) - sst2.collapsed("hour",iris.analysis.MIN)
  sst12_range = sst12.collapsed("hour",iris.analysis.MAX) - sst12.collapsed("hour",iris.analysis.MIN)
  amp12,peak12,mean12 = diurnal_amp_peak(P12)
  amp2,peak2,mean2 = diurnal_amp_peak(P2)
  peak_2  = peak2.copy()
  peak_2.data =np.ma.masked_array(peak_2.data,(amp2.data<1)+(mean2.data<2.5))
  peak_12  = peak12.copy()
  peak_12.data =np.ma.masked_array(peak_12.data,(amp12.data<1)+(mean12.data<2.5))
#
  diff = peak_2 - peak_12.regrid(peak_2,iris.analysis.Linear())
  diff = diff.data%24
  diff[diff>12] = diff[diff>12] - 24
  diff = peak_2.copy(data=diff)
#
  diff_sst= sst2_range-sst12_range.regrid(sst2_range,iris.analysis.Nearest())
  cmap1h@ = ListedColormap(sns.color_palette("turbo",12))
  cmap2 = ListedColormap(sns.color_palette('twilight_shifted',12))
  cmap3=ListedColormap(sns.color_palette("viridis",8))
  cmap4=ListedColormap(sns.color_palette("bwr",9))
  fig=plt.figure(figsize=(10,5))
  
  ax,p = [],[]
  for i,(cube,cmap,vmin,vmax,title,extend) in enumerate([(peak_2    ,cmap1,  0,24  ,"MC2 Hour of Max Precipitation",None),
                                     (peak_12    ,cmap1,  0,24,  "MC12 Hour of Max Precipitation",None),
                                     (diff       ,cmap4,-6,6,  "MC12 - MC12: Hour of Max Precipitation","Both"),
                                     (sst2_range ,cmap3,  0,0.8, "MC2 Diurnal SST Range",None),
                                     (sst12_range,cmap3,  0,0.8,  "MC12 Diurnal SST Range",None),
                                     (diff_sst   ,cmap4,-0.225,0.225,"MC2 - MC12 Diurnal SST Range",None)]):
    ax.append(plt.subplot(2,3,1+i,projection=ccrs.PlateCarree()))
    ax[-1].coastlines()
    plt.xlim(90,155)
    plt.ylim(-15,15)
    p.append(iplt.pcolormesh(cube,vmin=vmin,vmax=vmax,cmap=cmap))
    plt.title(title)
#
  mask = diff.copy(data=np.ma.masked_array(np.ones(diff.shape)))
  mask.data.mask = (diff.data <6)*(diff.data>-6)
  mask.data.mask += peak_12.data.mask
  iplt.pcolormesh(mask,axes=ax[2],cmap="binary",vmin=0,vmax=2)

  fig.colorbar(p[0],ax=ax[0],ticks=np.arange(0,25,4),orientation="horizontal")
  fig.colorbar(p[1],ax=ax[1],ticks=np.arange(0,25,4),orientation="horizontal")
  fig.colorbar(p[3],ax=ax[3],ticks=np.arange(0,.81,0.2),orientation="horizontal")
  fig.colorbar(p[4],ax=ax[4],ticks=np.arange(0,.81,0.2),orientation="horizontal")
  fig.colorbar(p[5],ax=ax[5],ticks=np.arange(-0.2,0.21,0.1),orientation="horizontal")
  c=fig.colorbar(p[2],ax=ax[2],ticks=np.arange(-6,7,3),orientation="horizontal")
  c.set_ticklabels(["Late"]+list(range(-8,9,4))+["Early"])
  fig.subplots_adjust(top=0.964,bottom=0.048,left=0.026,right=0.967,hspace=0.125,wspace=0.136)
  fig.savefig(figname)
  plt.show()
  import pdb;pdb.set_trace()


if __name__ == "__main__":
  plot([2015,2016],"/home/users/emmah/eval_figs/diurnal_cycle.png")



