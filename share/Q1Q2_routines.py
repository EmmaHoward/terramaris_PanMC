import matplotlib.pyplot as plt
from iris.aux_factory import HybridHeightFactory
import stratify
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np
import datetime as dt
from iris.analysis import calculus_centred,calculus
import os
import sys
from iris.coord_categorisation import add_hour,add_day_of_year

def interp_areaweighted(data,template):
# regrid all cubes in data to template grid
# bounds of 'template' are re-computed
  template.coord('latitude').bounds = None
  template.coord('longitude').bounds = None
  template.coord('longitude').guess_bounds()
  template.coord('latitude').guess_bounds()
  new  = iris.cube.CubeList([])
  for cube in data:
    cube.coord('longitude').coord_system = template.coord('longitude').coord_system
    cube.coord('latitude').coord_system = template.coord('latitude').coord_system
    try:
      cube.coord("longitude").guess_bounds()
      cube.coord("latitude").guess_bounds()
    except ValueError:
      1
    # regrid conservatively. cells that are at least 80% unmasked (due to orography) are retailed
    new.append(cube.regrid(template,iris.analysis.AreaWeighted(mdtol=0.2)))
  return new


def read_precip_MC108(dates,i,path,MC):
  alpha_w = iris.coords.AuxCoord(1,units='mm*m**2*kg**-1')
  year = dates[0].year - int(dates[0].month <6)
  ct = iris.Constraint(time=lambda t: dt.datetime(t.point.year,t.point.month,t.point.day) in dates)
  pp = iris.load("/gws/nopw/j04/terramaris/emmah/Q1Q2_analysis/%04d%02d_%s_coarsened_total_precip_MC%d.nc"%(year,(year+1)%100,MC,12*i),ct)[0]
  pp=pp*alpha_w
  pp.convert_units("mm")
  pp.rename("P")
  return pp

def read_Q1Q2_MC108(dates,i,path,MC):
  name = {"MC2":"2km","MC12":"N1280"}[MC]
  data = iris.load([path+"budgets/tma_%s_KPPcoupled_budget_%03d_%04d%02d%02d.nc"%(name,i*12,t.year,t.month,t.day) for t in dates])#.concatenate()
  for cube in data:
    cube.data = cube.data.astype("float")
  data=data.concatenate()
  Q1_e = data.extract("eddy_vertical_heat_flux_divergence")[0]
  Q2_e =-1000*data.extract("eddy_vertical_moisture_flux_divergence")[0]
  Q1_s = data.extract("subgrid_temperature_increment")[0]
  Q2_s =-1000*data.extract("subgrid_specific_humidity_increment")[0]
  Q1_s.data.mask = Q1_e.data.mask
  Q2_s.data.mask = Q2_e.data.mask
  Q1 = Q1_s - Q1_e
  Q2 = Q2_s - Q2_e
  #rho = data.extract("air_density")[0]
  Q1.rename("Q1")
  Q2.rename("Q2")
  return iris.cube.CubeList([Q1,Q2])#,rho])


def apply_mask(cubes,mask,coords=True):
  it,iy,ix = np.where(mask)
  new = {}
  for cube in cubes:
    if cube.ndim==4:
      assert cube[:,0].shape==mask.shape
      new[cube.name()]=cube.data[it,:,iy,ix]
    if cube.ndim==3:
      assert cube.shape==mask.shape
      new[cube.name()]=cube.data[it,iy,ix]
  if coords:
    nt,ny,nx = mask.shape
    tt,yy,xx = np.meshgrid(range(nt),range(ny),range(nx),indexing="ij")
    t = tt[it,iy,ix]
    y = yy[it,iy,ix]
    x = xx[it,iy,ix]
    return (new,t,y,x)
  else:
    return new


