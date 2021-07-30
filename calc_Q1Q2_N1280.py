#!/apps/contrib/jaspy/miniconda_envs/jaspy3.7/m3-4.5.11/envs/jaspy3.7-m3-4.5.11-r20181219/bin/python
#SBATCH -p short-serial
#SBATCH --array=[1-6]
#SBATCH -o /home/users/emmah/log/budget_15-%a.o
#SBATCH -e /home/users/emmah/log/budget_15-%a.e 
#SBATCH -t 01:00:00
#SBATCH --mem=24000



from scipy.fftpack import dct,idct
import stratify
import iris
from iris.experimental import equalise_cubes
import numpy as np
import datetime as dt
from iris.aux_factory import HybridHeightFactory
from iris.analysis import calculus_centred,calculus
import os


path_template = "/gws/nopw/j04/terramaris/emmah/um/u-bs742/processed/MC_15_orog.nc"
#path = "/gws/nopw/j04/terramaris/emmah/testrun/u-bs742/processed/terramaris_15km_MC_GA7/"

path = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/production_runs/201516_u-cf309/"

cx = iris.Constraint(longitude=lambda lon: 90<=lon<=155)
cy = iris.Constraint(latitude=lambda lat: -15<=lat<=15)

cx = iris.Constraint(longitude = lambda x: 90.5 < x < 154.4)
cy = iris.Constraint(latitude  =  lambda y: -15 < y < 15)

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
  q = budget.extract("specific_humidity")[0]
  mp1 = q/(1+-1*q)+1                        
  mp1.rename("m+1")

# calculate grid-spacing
  z_theta = dT_p[0].copy(data= dT_p.coord("altitude").points)
  z_theta.units="m"
  dz_theta = calculus.differentiate(z_theta,"model_level_number")
  dz_theta = q[0].copy(data=np.concatenate([[dT_p[:,0].coord("altitude").points-z_theta.coord("surface_altitude").points],dz_theta.data],axis=0))
  dz_theta.rename("theta_grid_spacing")
  dz_theta.units="m"

# adjust to moist density and remove 1e-5dz scaling from covariance terms. Ensure units are correct.
  wqrho = budget.extract("w_q_rho_dry_1e-5")[0]*1e5*mp1/dz_theta
  wqrho.units = "kg/m2/s"
  wqrho.rename("rho q w")
  wTrho = budget.extract("w_T_rho_dry_1e-5")[0]*1e5*mp1/dz_theta
  wTrho.units = "kg.K/m2/s"
  wTrho.rename("rho T w")

  data = iris.cube.CubeList([dT_p,dq_p,wqrho,wTrho]+budget.extract(["specific_humidity","air_temperature","upward_air_velocity"]))
  return data

def exner_rho(date,mp1):
# read rho from pa stream.
# if calculating potential temperature increment, uncomment code to read pressure and calculate exner pressure here
  rho = iris.load(path+"pa/tma_N1280_KPPcoupled_pa_density_%d%02d%02d.nc"%(date.year,date.month,date.day),cx&cy)[0]
#  p = iris.load(path+"pa/tma_N1280_KPPcoupled_pa_pressure_rho_grid_%d%02d%02d.nc"%(date.year,date.month,date.day),cx&cy)[0]
  for cube in [rho]:
    cube.remove_coord("forecast_period")
    cube.remove_coord("time")
    cube.add_dim_coord(mp1.coord("time"),0)
  rho = rho*mp1
#  p0 = iris.coords.AuxCoord(100000, long_name='P0', units='Pa')
#  pi = (p / p0)**(287.05/1005)
#  pi.rename("exner_pressure")
  rho.rename("air_density")
  data = iris.cube.CubeList([rho])
  return(data)
 
def regrid_z_asl(data):
  # convert from model level to above sea level data. data below orography is masked
  #z_t = data.extract("ddz_air_temperature")[0].coord('atmosphere_hybrid_height_coordinate').copy() # points to interpolate to: model t grid over the ocean
  z_t = np.array([2.0000000e+01, 5.3335999e+01, 1.0000000e+02,
                   1.6000000e+02, 2.3333600e+02, 3.2000000e+02,
                   4.2000000e+02, 5.3333600e+02, 6.6000000e+02,
                   8.0000000e+02, 9.5333600e+02, 1.1200000e+03,
                   1.3000000e+03, 1.4933361e+03, 1.7000000e+03,
                   1.9200000e+03, 2.1533359e+03, 2.4000000e+03,
                   2.6600000e+03, 2.9333359e+03, 3.2200000e+03,
                   3.5200000e+03, 3.8333359e+03, 4.1600000e+03,
                   4.5000000e+03, 4.8533359e+03, 5.2200000e+03,
                   5.6000000e+03, 5.9933359e+03, 6.4000000e+03,
                   6.8200000e+03, 7.2533442e+03, 7.7000400e+03,
                   8.1601362e+03, 8.6337041e+03, 9.1209043e+03,
                   9.6219600e+03, 1.0137232e+04, 1.0667248e+04,
                   1.1212736e+04, 1.1774704e+04, 1.2354504e+04,
                   1.2953904e+04, 1.3575160e+04, 1.4221144e+04,
                   1.4895432e+04, 1.5602456e+04, 1.6347608e+04,
                   1.7137424e+04, 1.7979729e+04, 1.8883840e+04,
                   1.9860775e+04]) 
#  z_t = np.arange(0,20000,100)
  new = iris.cube.CubeList()
  z_t = iris.coords.DimCoord(z_t,units="m", standard_name="altitude")
  for cube in data:
    z = cube.coord("altitude").points
    newdata = stratify.interpolate(z_t.points,z,cube.data.filled(np.nan),axis=1)
    new.append(cube[:,:len(z_t.points)].copy(data=newdata))
    for coord in ["model_level_number","sigma",'atmosphere_hybrid_height_coordinate']: 
      new[-1].remove_coord(coord)
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
      continue
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
  out = dcube/dz   out.rename("ddz_"+cube.name())
  out.rename("ddz_"+cube.name())
  return out


def main(date,i):
# calculate Q1 and Q2 for datetime date, coarsened over i*i N1280 grid-cells

  # read pc budget terms
  data = load_budget(date)            

  # calculate water vapour adjustment to density 
  q = data.extract("specific_humidity")[0]
  mp1 = q/(1+-1*q)+1                         
  mp1.rename("m+1")

  # read pa instantaneous terms
  data2 = exner_rho(date,mp1)   
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
  iris.save([dT_sub,dq_sub,dT_sg,dq_sg],path+"budgets/tma_N1280_KPPcoupled_budget_%03d_%04d%02d%02d.nc"%(i*12,date.year,date.month,date.day),zlib=True)


job =  int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
main(dt.datetime(2015,12,1)+dt.timedelta(job),9)

