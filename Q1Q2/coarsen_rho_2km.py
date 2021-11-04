#!/apps/contrib/jaspy/miniconda_envs/jaspy3.7/m3-4.5.11/envs/jaspy3.7-m3-4.5.11-r20181219/bin/python
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH --array=[7-17]
#SBATCH -o /home/users/emmah/log/budget_N1280/rho_2-%a.o
#SBATCH -e /home/users/emmah/log/budget_N1280/rho_2-%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=128000

from iris.aux_factory import HybridHeightFactory
from scipy.fftpack import dct,idct
import stratify
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np
import datetime as dt
from iris.aux_factory import HybridHeightFactory
from iris.analysis import calculus_centred,calculus
import os
import sys
from area_weighted_regrid_constant_altitude import AreaWeightedAltitude

t1 = dt.datetime(2015,12,13)

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1


#path = "/gws/nopw/j04/terramaris/emmah/testrun/u-bs742/processed/terramaris_15km_MC_GA7/"

def add_altitude_coord(cube):
        orog=iris.load(orogpath)[0]#.extract(cx&cy)
        orog_c = iris.coords.AuxCoord(orog.data,standard_name='surface_altitude',units='m')
        cube.add_aux_coord(orog_c,[2,3])
        a = HybridHeightFactory(cube.coord('level_height'),\
                          cube.coord('sigma'),\
                          cube.coord('surface_altitude'))
        a.rename("altitude")
        cube.add_aux_factory(a)


path = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/%04d%02d_u-cc339/"%(year,(year+1)%100)
orogpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/orog.pp"

cx = iris.Constraint(longitude=lambda lon: 90<=lon<=155)
cy = iris.Constraint(latitude=lambda lat: -15<=lat<=15)


def Tgrid_to_rhogrid(cube,target,lwr_bdy):
  if cube.ndim==4:
    data=np.concatenate([lwr_bdy.data.reshape(cube[:,:1].shape),cube.data],axis=1)
    new = target.copy(data = (data[:,1:]+data[:,:-1])/2)
    new.units=cube.units
    new.rename(cube.name())
    return new
  if cube.ndim==3:
    data=np.concatenate([lwr_bdy.data.reshape(cube[:1].shape),cube.data],axis=1)
    new = target.copy(data = (data[1:]+data[:-1])/2)
    new.units=cube.units
    new.rename(cube.name())
    return new



def interp2km(cube):
  path_template = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"
  template = iris.load(path_template)[0][54:-55,42:-42] # bounds extract inner domain
#  orog_c = iris.coords.AuxCoord(template.data,long_name='surface_altitude',units='m')
#  template.add_aux_coord(orog_c,(0,1))
  if not template.coord('longitude').has_bounds():
    template.coord('longitude').guess_bounds()
    template.coord('latitude').guess_bounds()
  if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()
        cube.coord('latitude').guess_bounds()
  #template = template.extract(cx&cy)
  new = iris.cube.CubeList()
  orog_c = cube[0,0].copy(data = cube.coord("surface_altitude").points).regrid(template,iris.analysis.AreaWeighted())
  orog_c = iris.coords.AuxCoord(orog_c.data,long_name='surface_altitude',units='m')
  template.add_aux_coord(orog_c,(0,1))
  if cube.data.dtype == np.double:
           cube.data = cube.data.astype(float)
  if cube.ndim==4:
    for t in range(cube.shape[0]):
      print(cube[t])
      new.append(cube[t].regrid(template,AreaWeightedAltitude("altitude","surface_altitude")))
  new = new.merge_cube()
  return new

def exner_rho(date):
# read rho from pa stream.
  rhod = iris.load(path+"pa/tma_2km_KPPcoupled_pa_density_%d%02d%02d_%02d.nc"%(date.year,date.month,date.day,date.hour),"air_density")[0]
  q = iris.load(path+"pa/tma_2km_KPPcoupled_pa_specific_humidity_%d%02d%02d_%02d.nc"%(date.year,date.month,date.day,date.hour),"specific_humidity")[0]
  ct = iris.Constraint(time=lambda t: t.point.hour in [tt.hour for tt in q.coord("time").units.num2date(q.coord("time").points)])
  qs = iris.load(path+"pb/tma_2km_KPPcoupled_pb_surf_inst_%d%02d%02d.nc"%(date.year,date.month,date.day),"specific_humidity")[0]#.extract(ct)
  qs = qs.extract(ct)
  q_r = Tgrid_to_rhogrid(q,rhod,qs)
  mp1 = q_r/(1+-1*q_r)+1  # convert dry to moist density                      
  rho = rhod*mp1
  rho.rename(rhod.name())
  add_altitude_coord(rho)
  rho = interp2km(rho)
  iris.save(rho,"/work/scratch-nopw/emmah/rho/rho_%d%02d%02d_%02d.nc"%(date.year,date.month,date.day,date.hour))

job =int(os.environ["SLURM_ARRAY_TASK_ID"]) -1 
for i in range(0,24,4):
  exner_rho(t1+dt.timedelta(job+i/24))

