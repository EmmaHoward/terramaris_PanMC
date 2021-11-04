#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python

import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call
from iris.aux_factory import HybridHeightFactory
import sys


#cx = iris.Constraint(longitude=lambda x: 100<=x<=102)
#cy = iris.Constraint(latitude=lambda x: -5<=x<=-3)
#cz = iris.Constraint(model_level_number = lambda x: x<10)

# Constraints on cell methods (ie UM time profiles)
c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="mean"))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not cube.cell_methods)

# system paths
orogpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/orog.pp"
outpath = "/work/scratch-nopw/emmah/postproc_2km_fc/"
finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/"
# re-initialisation timesteps
#if t1.day in [1,2]:
#  reinit_step_pp = {"pa":4,"pb":4,"pc":4,"pd":6} # input, in hours
#else:
reinit_step_pp = {"pa":1,"pb":4,"pc":1,"pd":6} # input, in hours
reinit_step_nc = {"pa":1,"pb":6,"pc":1,"pd":6}   # output, in days
reinit_step_nc = {"pa":4,"pb":24,"pc":4,"pd":24}   # output, in days

reinit_step_pa = {"density":4,"potential_temperature":4,"pressure_rho_grid":4,"pressure_theta_grid":4,\
                  "upward_air_velocity":4,"eastward_wind":4,"northward_wind":4,"specific_humidity":4,\
                  "radar_reflectivity":24,"cloud_fraction":24,"qcf":24,"qcl":24,"rain":24,"graupel":24}

from file_split_2km import file_variables_full as file_variables

# variable names to be corrected
names = {
  "x_wind":"eastward_wind",  # note: terramaris grids are not rotated
  "y_wind":"northward_wind",
# pa
  "m01s05i250":"atmosphere_updraft_convective_mass_flux",
  "m01s05i251":"atmosphere_downdraft_convective_mass_flux",
# pb
  "m01s01i202":"surface_net_downward_shortwave_flux",
  "m01s03i241":"surface_moisture_flux",
  "m01s03i319":"canopy_height",
  "m01s03i359":"height_of_diagnostic_parcel_top",
  "m01s03i462":"stomatal_conductance",
  "m01s03i539":"transpiration_amount",
  "m01s03i465":"explicit_friction_velocity",
  "m01s03i476":"combined_boundary_layer_type",
  "m01s09i202":"very_low_type_cloud_area_fraction",
  "m01s30i403":"total_column_dry_mass",
  "m01s30i461":"atmosphere_mass_content_of_water_vapor",
  "m01s30i462":"Vertical_integral_of_eastward_water_vapour_flux",
  "m01s30i463":"Vertical_integral_of_northward_water_vapour_flux",
# pc
  "m01s01i161":"change_over_time_in_air_temperature_due_to_shortwave_radiation_noPC2",
  "m01s02i161":"change_over_time_in_air_temperature_due_to_longwave_radiation_noPC2",
  "m01s30i034":"w_T_rho_1e-5",
  "m01s30i035":"w_q_rho_1e-5",
# pd
  "m01s15i229":"potential_vorticity",
  "Heavyside function on pressure levels":"Heaviside function on pressure levels"
  }

# units to be corrected
units = {
# pa
  "m01s05i250":"Pa/s",
  "m01s05i251":"Pa/s",
# pb 
  "m01s01i202":"W/m2",
  "m01s03i319":"m",
  "m01s03i359":"m",
  "m01s03i462":"m/s",
  "m01s03i465":"m/s",
  "m01s03i476":"1",
  "m01s03i539":"kg/m2",
  "m01s09i202":"1",
  "m01s30i403":"kg/m2",
  "m01s30i461":"kg/m2",
  "m01s30i462":"kg/m/s",
  "m01s30i463":"kg/m/s",
# pc
  "m01s01i161":"K/hr",
  "m01s02i161":"K/hr",
  "m01s30i034":"kg.K/m2/s",#only correct after dividing by dz
  "m01s30i035":"kg/m2/s",  #only correct after dividing by dz
  "change_over_time_in_air_temperature_due_to_advection":"K/hr",
  "change_over_time_in_specific_humidity_due_to_advection":"1/hr",
  "change_over_time_in_air_temperature":"K/hr",
  "change_over_time_in_specific_humidity":"1/hr",
# pd
  "m01s15i229":"m2.K/s/kg",
  "divergence_of_wind":"10e6 s-1",
  "atmosphere_relative_vorticity":"10e6 s-1"
}

from area_weighted_regrid_constant_altitude import AreaWeightedAltitude
def interp(cube):
  path_template = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"
  template = iris.load(path_template)[0][54:-55,42:-42]#.extract(cx&cy) # bounds extract inner domain
#  orog_c = iris.coords.AuxCoord(template.data,standard_name='surface_altitude',units='m')
#  template.add_aux_coord(orog_c,(0,1))
  if not template.coord('longitude').has_bounds():
    template.coord('longitude').guess_bounds()
    template.coord('latitude').guess_bounds()
  if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()
        cube.coord('latitude').guess_bounds()
#  template=template.extract(cx&cy)
  new = iris.cube.CubeList()
  if cube.ndim==4:
    orog_c = cube[0,0].copy(data = cube.coord("surface_altitude").points).regrid(template,iris.analysis.AreaWeighted())
    orog_c = iris.coords.AuxCoord(orog_c.data,long_name='surface_altitude',units='m')
    template.add_aux_coord(orog_c,(0,1))
    if cube.data.dtype == np.double:
           cube.data = cube.data.astype(float)
    for t in range(cube.shape[0]):
      print(cube[t])
      new.append(cube[t].regrid(template,AreaWeightedAltitude("altitude","surface_altitude")))
  if cube.ndim==3:
    orog_c = cube[0].copy(data = cube.coord("surface_altitude").points).regrid(template,iris.analysis.AreaWeighted())
    orog_c = iris.coords.AuxCoord(orog_c.data,long_name='surface_altitude',units='m')
    template.add_aux_coord(orog_c,(0,1))
    if cube.data.dtype == np.double:
           cube.data = cube.data.astype(float)
    new.append(cube.regrid(template,AreaWeightedAltitude("altitude","surface_altitude")))
  new = new.merge_cube()
  return new


def add_altitude_coord(cube):
        orog=iris.load(orogpath)[0]#.extract(cx&cy)
        orog_c = iris.coords.AuxCoord(orog.data,standard_name='surface_altitude',units='m')
        if cube.ndim==4:
          cube.add_aux_coord(orog_c,[2,3])
        if cube.ndim==3:
          cube.add_aux_coord(orog_c,[1,2])
        a = HybridHeightFactory(cube.coord('level_height'),\
                          cube.coord('sigma'),\
                          cube.coord('surface_altitude'))
        a.rename("altitude")
        cube.add_aux_factory(a)

def get_wxrho_factor(data):
        dthetadt = data.extract("change_over_time_in_air_temperature")[0] # data on theta grid
        q = data.extract("specific_humidity")[0]
        mp1 = q/(1+-1*q)+1  # convert dry to moist density                      
        mp1.rename("m+1")
        add_altitude_coord(dthetadt)
        z_T = dthetadt.coord("altitude").points
        dz_T = z_T[1:]-z_T[:-1]
        dz_T = np.concatenate([[z_T[0]],  dz_T]) 
        fac = mp1.data/dz_T
        return fac,dz_T,mp1

import matplotlib.pyplot as plt
import iris.plot as iplt 
def main():
  data=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/2015/20151207T0000Z/tm2a_pc20151212_00")
  rho=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/2015/20151207T0000Z/tm2a_pa20151212_00","air_density")
  for cube in data:
    if cube.name() in units.keys():
        # correct units
        cube.units = units[cube.name()]
    if cube.name() in names.keys():
        # correct variable names
        cube.rename(names[cube.name()])
  fac,dz_T,mp1 = get_wxrho_factor(data)  
  add_altitude_coord(rho[0])
  earth_radius = 6371229.0
  z_r = rho[0].coord("altitude").points
  rho[0].remove_aux_factory(rho[0].aux_factories[0])
  import pdb;pdb.set_trace()
  wtrho_tot = data.extract("w_T_rho_1e-5")[0][:,200,200]*1e5/dz_T[:,200,200]*mp1[:,200,200]*(earth_radius/(z_r[:,200,200]+earth_radius))**2
  wtrho_mean = data.extract("air_temperature")[0][:,200,200]*rho[0][:,200,200]*data.extract("upward_air_velocity")[0][:,200,200]*mp1[:,200,200]
  wqrho_tot = data.extract("w_q_rho_1e-5")[0][:,200,200]*1e5/dz_T[:,200,200]*mp1[:,200,200]*(earth_radius/(z_r[:,200,200]+earth_radius))**2
  wqrho_mean = data.extract("specific_humidity")[0][:,200,200]*rho[0][:,200,200]*data.extract("upward_air_velocity")[0][:,200,200]*mp1[:,200,200]

main()
