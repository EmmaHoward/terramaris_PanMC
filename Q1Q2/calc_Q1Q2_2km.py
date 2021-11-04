#!/apps/contrib/jaspy/miniconda_envs/jaspy3.7/m3-4.5.11/envs/jaspy3.7-m3-4.5.11-r20181219/bin/python
#  SBATCH -p test
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH --array=[1-36]
#SBATCH -o /home/users/emmah/log/budget_N1280/budget_15-%a.o
#SBATCH -e /home/users/emmah/log/budget_N1280/budget_15-%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=48000



import stratify
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np
import datetime as dt
from iris.aux_factory import HybridHeightFactory
from iris.analysis import calculus_centred,calculus
import os
import sys

#t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))
t1 = dt.datetime(2015,12,1)

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1




#path = "/gws/nopw/j04/terramaris/emmah/testrun/u-bs742/processed/terramaris_15km_MC_GA7/"

path = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/%04d%02d_u-cc339/"%(year,(year+1)%100)
#path = "/work/scratch-nopw/emmah/postproc_2km_fc/"
cx = iris.Constraint(longitude=lambda lon: 90<=lon<=155)
cy = iris.Constraint(latitude=lambda lat: -15<=lat<=15)

#cx = iris.Constraint(longitude = lambda x: 90.5 < x < 154.4)
#cy = iris.Constraint(latitude  =  lambda y: -15 < y < 15)

#cx = iris.Constraint(longitude=lambda lon: 120<=lon<=124)
#cy = iris.Constraint(latitude=lambda lat: -4<=lat<=0)


coord_names ={"height based hybrid coeffient a":'level_height',
              "height based hybrid coeffient b":"sigma"}




def load_budget(date):
# read in all terms from pc stream
  budget = iris.load(path+"pc/*%d%02d%02d.nc"%(date.year,date.month,date.day),cx&cy)
  for cube in budget:
    try:
      cube.remove_coord("forecast_period")
    except: 
      continue
# calculate subgrid-scale heating and drying terms as residuals from total, advection and radiative terms
  dT_p =  budget.extract("change_over_time_in_air_temperature")[0] \
         - budget.extract("change_over_time_in_air_temperature_due_to_advection")[0] \
         - budget.extract("change_over_time_in_air_temperature_due_to_longwave_radiation_noPC2")[0] \
         - budget.extract("change_over_time_in_air_temperature_due_to_shortwave_radiation_noPC2")[0]
  dT_p.rename("subgrid_temperature_increment")
  dq_p =  budget.extract("change_over_time_in_specific_humidity")[0] \
         - budget.extract("change_over_time_in_specific_humidity_due_to_advection")[0]
  dq_p.rename("subgrid_specific_humidity_increment")

# calculate water vapour adjustment to density
#  q = budget.extract("specific_humidity")[0]
#  mp1 = q/(1+-1*q)+1                        
#  mp1.rename("m+1")

# remove 1e-5dz scaling from covariance terms. Ensure units are correct.
  wqrho = budget.extract("w_q_rho_1e-5")[0]*1e5
  print(wqrho.units)
  print("kg/m2/s")
#  wqrho.units = "kg/m2/s"
  wqrho.rename("rho q w")
  wTrho = budget.extract("w_T_rho_1e-5")[0]*1e5
#  wTrho.units = "kg.K/m2/s"
  wTrho.rename("rho T w")

  data = iris.cube.CubeList([dT_p,dq_p,wqrho,wTrho]+budget.extract(["specific_humidity","air_temperature","upward_air_velocity"]))
  #for cube in data:
  #  for coord in ["longitude","latitude","surface_altitude"]:
  #     cube.coord(coord).var_name=None
  return data

def exner_rho(date):
# read rho from pa stream.
# if calculating potential temperature increment, uncomment code to read pressure and calculate exner pressure here

  rho = iris.load("/work/scratch-nopw/emmah/rho/rho_%d%02d%02d_*.nc"%(date.year,date.month,date.day),cx&cy)
  if not (date.day==1 and date.month==12):
    datem1 = date+dt.timedelta(-1)
    rho += iris.load("/work/scratch-nopw/emmah/rho/rho_%d%02d%02d_20.nc"%(datem1.year,datem1.month,datem1.day),cx&cy)
  equalise_attributes(rho)
  rho = rho.extract("air_density").concatenate_cube()
  a = HybridHeightFactory(rho.coord('atmosphere_hybrid_height_coordinate'),\
                          rho.coord('sigma'),\
                          rho.coord('surface_altitude'))
  rho.add_aux_factory(a)
  if date.day==1 and date.month==12:
    rho0 = rho[0]
    rho0.coord("time").points = 720.5
    rho0.coord("time").bounds = [[720,721]]
    rho_roll = rho.rolling_window("time",iris.analysis.MEAN,2)
    rho_roll.data=rho_roll.data.astype("float32")
    rho0.cell_methods = rho_roll.cell_methods
    rho_new = iris.cube.CubeList([rho0]+[rho_roll[i] for i in range(23)]).merge_cube()
    rho = rho_new
  else:
    rho = rho[3:]
    rho = rho.rolling_window("time",iris.analysis.MEAN,2)
#  p = iris.load(path+"pa/tma_N1280_KPPcoupled_pa_pressure_rho_grid_%d%02d%02d.nc"%(date.year,date.month,date.day),cx&cy)[0]
#  for cube in [rho]:
#    cube.remove_coord("forecast_period")
#    cube.remove_coord("time")
#    cube.add_dim_coord(mp1.coord("time"),0)
#  rho = rho*mp1
#  p0 = iris.coords.AuxCoord(100000, long_name='P0', units='Pa')
#  pi = (p / p0)**(287.05/1005)
#  pi.rename("exner_pressure")
  rho.rename("air_density")
  data = iris.cube.CubeList([rho])
  return(data)
 
def regrid_z_asl(data):
  # convert from model level to above sea level data. data below orography is masked
  z_t = data.extract("subgrid_temperature_increment")[0].coord('atmosphere_hybrid_height_coordinate').points.copy() # points to interpolate to: model t grid over the ocean
  new = iris.cube.CubeList()
  z_t = iris.coords.DimCoord(z_t,units="m", standard_name="altitude")
  for cube in data:
    a = HybridHeightFactory(cube.coord('atmosphere_hybrid_height_coordinate'),\
                          cube.coord('sigma'),\
                          cube.coord('surface_altitude'))
    a.rename("altitude")
    #cube.add_aux_factory(a)
    z = a.make_coord(cube.coord_dims).points
    newdata = stratify.interpolate(z_t.points,z,cube.data.filled(np.nan),axis=1)
    new.append(cube[:,:len(z_t.points)].copy(data=newdata))
    for coord in ["model_level_number","sigma",'atmosphere_hybrid_height_coordinate',"surface_altitude"]: 
      if coord in [c.name() for c in new[-1].coords()]:
        new[-1].remove_coord(coord)
    if len(new[-1].aux_factories) > 0:
      new[-1].remove_aux_factory(new[-1].aux_factories[0])
    new[-1].add_dim_coord(z_t,1)
  return new

def interp(data,template):
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


def ddz(cube):
  # calculate vertical derivative using first order centred difference
  out = calculus_centred.differentiate(cube,"altitude")
  # return to original grid (adds a masked point above and below new data)
  out = out.interpolate([("altitude",cube.coord("altitude").points)],iris.analysis.Linear(extrapolation_mode='mask'))
  out.rename("ddz_"+cube.name())
  return out

def ddz_(cube):
  # legacy: calculate vertical derivative from a sigma-coordinate grid
  z = cube[0].copy(data=cube.coord("altitude").points)
  z.units="m"
  dz = calculus.differentiate(z,"model_level_number")
  dcube = calculus.differentiate(cube,"model_level_number")
  out = dcube/dz   
  out.rename("ddz_"+cube.name())
  return out


def main(date,i):
# calculate Q1 and Q2 for datetime date, coarsened over i*i N1280 grid-cells

  # read pc budget terms
  data = load_budget(date)            

  # read pa instantaneous terms
  data2 = exner_rho(date)   
  data=data+data2
  print(data)
  # regrid vertically from sigma levels to altitude above sea-level. New grid is theta-levels above ocean points
  data=regrid_z_asl(data)

  # regrid horizontally to template grid centred on every ith grid-cell
  template = data[0][0,0,i//2:-i//2-1:i,i//2:-i//2-1:i]
  data = interp(data,template)

  # compute reynolds averaged large-scale fluxes
  rho_w_q_means = data.extract("specific_humidity")[0]*data.extract("upward_air_velocity")[0]*data.extract("air_density")[0]
  rho_w_q_means.rename("rho q w means")
  data.append(rho_w_q_means)
  rho_w_T_means = data.extract("air_temperature")[0]*data.extract("upward_air_velocity")[0]*data.extract("air_density")[0]
  rho_w_T_means.rename("rho T w means")
  data.append(rho_w_T_means)
  # compute vertical derivatives of fluxes
  for cube in data.extract(["rho T w", "rho q w","rho T w means","rho q w means"]):
    data.append(ddz(cube))

  # weight by density
  dT_tot = data.extract("ddz_rho T w")[0]/data.extract("air_density")[0]
  dq_tot = data.extract("ddz_rho q w")[0]/data.extract("air_density")[0]
  dT_tot.rename("total_vertical_heat_flux_divergence")
  dq_tot.rename("total_vertical_moisture_flux_divergence")

  dT_mean = data.extract("ddz_rho T w means")[0]/data.extract("air_density")[0]
  dq_mean = data.extract("ddz_rho q w means")[0]/data.extract("air_density")[0]
  dT_mean.rename("large-scale_vertical_heat_flux_divergence")
  dq_mean.rename("large-scale_vertical_moisture_flux_divergence")

  # convert units
  dT_mean.convert_units("K/hr")
  dT_tot.convert_units("K/hr")
  dq_mean.convert_units("1/hr")
  dq_tot.convert_units("1/hr")

  # compute reynolds averaged eddy terms
  dT_sub = dT_tot - dT_mean
  dq_sub = dq_tot - dq_mean
  dT_sub.rename("eddy_vertical_heat_flux_divergence")
  dq_sub.rename("eddy_vertical_moisture_flux_divergence")

  # extract subgrid-scale terms
  dT_sg = data.extract("subgrid_temperature_increment")[0]
  dq_sg = data.extract("subgrid_specific_humidity_increment")[0]
  dT_sg.units = "K/hr"
  dq_sg.units = "1/hr"
  # save to file
  #iris.save([dT_mean,dT_tot,dT_sub,dq_mean,dq_tot,dq_sub,dT_sg,dq_sg],path+"budgets/tma_2km_KPPcoupled_budget_%03d_%04d%02d%02d.nc"%(i*12,date.year,date.month,date.day),zlib=True)
  out = iris.cube.CubeList([dT_sub,dq_sub,dT_sg,dq_sg,data.extract("air_density")[0]])
  for cube in out:
    cube.data = cube.data.as_type(float)
  iris.save(out,path+"budgets/tma_2km_KPPcoupled_budget_%03d_%04d%02d%02d.nc"%(i*12,date.year,date.month,date.day),zlib=True)


job =  int(os.environ["SLURM_ARRAY_TASK_ID"]) -1 
main(t1+dt.timedelta(job),9)

