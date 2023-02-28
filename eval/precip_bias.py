#
#  Generates figure 4
#

import iris
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
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/precip/"
MC2_path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"
ref_path = "/gws/nopw/j04/klingaman/datasets/GPM_IMERG/monthly/"
cx = iris.Constraint(longitude = lambda x: 84 < x < 161)
cy = iris.Constraint(latitude  =  lambda y: -21 < y < 21)


def load(year,MC,path):
  """
  Loads precipitation data (diurnal means on MC12 grid)
  Parameters:
  =========
  year: start year of season to load 
  MC: MC2 or MC12
  path: input data path
  """
  if MC =="MC2":
    data=iris.load([path+"%04d%02d_diurnal_precip_coarsened.nc"%(year+int(month<6),month) for month in [12,1,2]])
  else:
    data=iris.load([path+"%04d%02d_diurnal_precip.nc"%(year+int(month<6),month) for month in [12,1,2]])
  equalise_attributes(data)
  data=data.concatenate()
  if MC=="MC12":
    P = (data[0]+data[1]+data[2]+data[3])*4*24
  else:
    P = (data[0])*4*24
  P.units="mm/day"
  try:
    P.coord("longitude").guess_bounds()
    P.coord("latitude").guess_bounds()
  except ValueError:
    1
  return P

def GPM(year,template):
  # loads GPM reference data
  # year: start year of season to load
  # template: grid to interpolate data onto
  gpm = iris.load([ref_path+"%04d/3B-MO.MS.MRG.3IMERG.%04d%02d01-S000000-E235959.%02d.V06B.nc"%(year+int(month<6),year+int(month<6),month,month) for month in [12,1,2]],cx&cy).concatenate_cube()
  gpm.coord("longitude").guess_bounds()
  gpm.coord("latitude").guess_bounds()
  gpm.coord("longitude").coord_system = template.coord("longitude").coord_system
  gpm.coord("latitude").coord_system = template.coord("latitude").coord_system
  gpm=gpm.regrid(template,iris.analysis.AreaWeighted())
  gpm.convert_units("mm/day")
  return gpm

def load_all(years12,years2):
  # load all data into memory 
  # years12: years to load for MC12 and the reference
  # years to load for MC2
  P2 = iris.cube.CubeList() #MC2
  P12 = iris.cube.CubeList() #MC12
  Pref = iris.cube.CubeList() #reference
  for year in years12:
    P12.append(load(year,"MC12",MC12_path))
    Pref.append(GPM(year,P12[-1]))
    P12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    Pref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  for year in years2:
    P2.append(load(year,"MC2",MC2_path))
    P2[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  equalise_attributes(P2)
  P12 = P12.concatenate_cube()
  P2 = P2.concatenate_cube()
  Pref = Pref.concatenate_cube()
  # regrid MC2 data to MC12
  P2 = P2.regrid(P12,iris.analysis.AreaWeighted(mdtol=0.5))
  return P12,P2,Pref


def main(years12,years2,figname=None):
  """
  Generate precipitation bias plots
  """
  # read data
  P12,P2,Pref = load_all(years12,years2)
  # time means
  P2= P2.collapsed("time",iris.analysis.MEAN)
  P12= P12.collapsed("time",iris.analysis.MEAN)
  Pref=Pref.collapsed("time",iris.analysis.MEAN)
  # area constraints for printing rainfall means away from boundaries
  cy2 = iris.Constraint(latitude=lambda x: -13<=x<=13)
  cx2 = iris.Constraint(longitude=lambda x: 95<=x<=150)
  # print rainfall means away from the boundaries for quoting in paper text
  print(P2.extract(cx2&cy2).collapsed(['latitude','longitude'],iris.analysis.MEAN).data)
  print(P12.extract(cx2&cy2).collapsed(['latitude','longitude'],iris.analysis.MEAN).data)
  print(Pref.extract(cx2&cy2).collapsed(['latitude','longitude'],iris.analysis.MEAN).data)
  if figname:
    # plot figure
    fig=bias_plots([P2,P12,Pref], ["MC2","MC12","GPM-IMERG"] ,cmap1="cmo.rain",cmap2="bwr_r",projection=ccrs.PlateCarree(),nmax=12,above=True,mark_inner=True,minmaxN1rangeN2=(0,25,11,11.25,12))
    fig.set_figwidth(12)
    fig.set_figheight(7)
    fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
    fig.suptitle("Precipitation (mm/day)")
    fig.tight_layout()
    fig.savefig(figname)

