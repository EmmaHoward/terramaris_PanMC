#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH --array=[15]
#SBATCH -o /home/users/emmah/log/postproc_2km/pp2_%a.o
#SBATCH -e /home/users/emmah/log/postproc_2km/pp2_%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=24000

import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call
from iris.aux_factory import HybridHeightFactory
year=2015
t1 = dt.datetime(2015,12,1)

region = {"pe":"Java","pf":"Bengkulu"}

# Constraints on cell methods (ie UM time profiles)
c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method  in ["mean","sum"]))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not (cube.cell_methods) or (cube.cell_methods[0].method=="point"))

# system paths
path = "/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/%04d/%04d%02d%02dT0000Z/"%(year,t1.year,t1.month,t1.day)
orogpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/orog.pp"
outpath = "/work/scratch-nopw/emmah/postproc_2km_fc/"
finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/"
# re-initialisation timesteps
reinit_step_pp = 1# {"pa":4,"pb":24,"pc":6,"pd":12} # input, in hours
reinit_step_nc = 1#{"pa":1,"pb":6,"pc":1,"pd":6}   # output, in days

# files and variables in the form:
# file name: (stream, [valid cell methods], [variable names])
file_variables = {
         "cloud_fraction"       :([c_inst],["cloud_area_fraction_in_atmosphere_layer",                               # 1
                                         "bulk_cloud_fraction_in_atmosphere_layer",
                                         "liquid_cloud_fraction_in_atmosphere_layer",
                                         "frozen_cloud_fraction_in_atmosphere_layer",
                                         "ice_aggregate_fraction_in_atmosphere_layer"]),

         "radar_reflectivity"   :([c_inst],["radar_reflectivity_due_to_graupel",                                     # 2
                                          "radar_reflectivity_due_to_cloud_ice",
                                          "radar_reflectivity_due_to_cloud_snow",
                                          "radar_reflectivity_due_to_rain",
                                          "radar_reflectivity_due_to_cloud_liquid"]),
         "density"              :([c_inst],["air_density"]),                                                         # 3
         "potential_temperature":([c_inst],["air_potential_temperature"]),                                           # 4
         "pressure_rho_grid"    :([c_inst],["m01s00i407"]),                                                          # 5
         "pressure_theta_grid"  :([c_inst],["m01s00i408"]),                                                          # 6
         "mixing_ratios"        :([c_inst],["mass_fraction_of_cloud_ice_in_air",                                     # 7
                                                 "mass_fraction_of_cloud_liquid_water_in_air",
                                                 "mass_fraction_of_graupel_in_air",
                                                 "mass_fraction_of_rain_in_air"]),
         "specific_humidity"    :([c_inst],["specific_humidity"]),                                                   # 8
         "winds"                :([c_inst],["eastward_wind",                                                         # 9
                                                 "northward_wind",
                                                 "upward_air_velocity"]), 
#
         "atmos":                ([c_mean,c_inst],["atmosphere_boundary_layer_thickness",                            #10
                                        "toa_incoming_shortwave_flux",
                                        "toa_outgoing_shortwave_flux",
                                        "toa_outgoing_shortwave_flux_assuming_clear_sky",
                                        "toa_outgoing_longwave_flux",
                                        "toa_outgoing_longwave_flux_assuming_clear_sky",
                                        "Turbulent mixing height after boundary layer",
                                        "height_of_diagnostic_parcel_top",
                                        "combined_boundary_layer_type"]),
         "integral_inst":        ([c_inst],["Stable boundary layer indicator",                                       #11
                                        "Stratocumulus over stable boundary layer indicator",
                                        "Well-mixed boundary layer indicator",
                                        "Decoupled stratocumulus not over cumulus indicator",
                                        "Decoupled stratocumulus over cumulus indicator",
                                        "Cumulus capped boundary layer indicator",
                                        "cloud_area_fraction_assuming_random_overlap",
                                        "cloud_area_fraction_assuming_maximum_random_overlap",
                                        "total_column_dry_mass",
                                        "atmosphere_mass_per_unit_area",
                                        "atmosphere_cloud_liquid_water_content",
                                        "atmosphere_cloud_ice_content",
                                        "atmosphere_mass_content_of_water_vapor"]),
         "integral_agg":         ([c_mean],[                                                                         #12
                                        "very_low_type_cloud_area_fraction",
                                        "low_type_cloud_area_fraction",
                                        "medium_type_cloud_area_fraction",
                                        "high_type_cloud_area_fraction",
                                        "Vertical_integral_of_eastward_water_vapour_flux",
                                        "Vertical_integral_of_northward_water_vapour_flux"]),
         "rainfall":                 ([c_mean,c_inst],["stratiform_rainfall_amount",                                 #13
                                        "stratiform_rainfall_flux"]),
         "surf_inst":            ([c_inst],["surface_temperature",                                                   #14
                                        "surface_air_pressure",
                                        "eastward_wind",
                                        "northward_wind",
                                        "air_temperature",
                                        "specific_humidity",
                                        "relative_humidity",
                                        "dew_point_temperature",
                                        "explicit_friction_velocity"]),
         "surf_agg":             ([c_mean,c_max,c_inst],[
                                        "m01s01i201",                                                                #15
                                        "m01s01i202",
                                        "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky",
                                        "surface_upwelling_shortwave_flux_in_air_assuming_clear_sky",
                                        "surface_downwelling_shortwave_flux_in_air",
                                        "surface_net_downward_longwave_flux",
                                        "surface_downwelling_longwave_flux_in_air",
                                        "surface_downwelling_longwave_flux_in_air_assuming_clear_sky",
                                        "surface_upward_sensible_heat_flux",
                                        "surface_upward_water_flux",
                                        "Evaporation flux from open sea",
                                        "surface_upward_latent_heat_flux",
                                        "surface_downward_eastward_stress",
                                        "surface_downward_northward_stress",
                                        "wind_speed_of_gust",
                                        "number_of_lightning_flashes"])}


names = {
# multiple
  "x_wind":"eastward_wind",
  "y_wind":"northward_wind",
  "NUMBER OF LIGHTNING FLASHES":"number_of_lightning_flashes",
  "ICE AGGREGATE FRACTION":"ice_aggregate_fraction_in_atmosphere_layer",
  "TOTAL RADAR REFLECTIVITY 3D":"radar_reflectivity_due_to_all_hydrometeor_species",
  "GRAUPEL RADAR REFLECTIVITY   (dBZ)":"radar_reflectivity_due_to_graupel",
  "ICE AGG RADAR REFLECTIVITY  (dBZ)":"radar_reflectivity_due_to_cloud_ice",
  "ICE CRY RADAR REFLECTIVITY (dBZ)":"radar_reflectivity_due_to_cloud_snow",
  "RAIN RADAR REFLECTIVITY (dBZ)":"radar_reflectivity_due_to_rain",
  "LIQ. CLOUD RADAR REFLECTIVITY (dBZ)":"radar_reflectivity_due_to_cloud_liquid",
  "BULK CLOUD FRACTION IN EACH LAYER":"bulk_cloud_fraction_in_atmosphere_layer",
  "LIQUID CLOUD FRACTION IN EACH LAYER":"liquid_cloud_fraction_in_atmosphere_layer",
  "FROZEN CLOUD FRACTION IN EACH LAYER":"frozen_cloud_fraction_in_atmosphere_layer",
  "GRAUPEL AFTER TIMESTEP":"mass_fraction_of_graupel_in_air",
#  "m01s05i250":"atmosphere_updraft_convective_mass_flux",
#  "m01s05i251":"atmosphere_downdraft_convective_mass_flux",
  "m01s01i202":"surface_net_downward_shortwave_flux",
  "m01s03i241":"surface_upward_water_flux",
  "m01s03i319":"canopy_height",
  "m01s03i359":"height_of_diagnostic_parcel_top",
  "m01s03i462":"stomatal_conductance",
  "m01s03i465":"explicit_friction_velocity",
  "m01s03i476":"combined_boundary_layer_type",
#  "surface_water_flux":"surface_upward_water_flux",
  "m01s09i202":"very_low_type_cloud_area_fraction",
  "m01s30i403":"total_column_dry_mass",
  "m01s30i461":"atmosphere_mass_content_of_water_vapor",
  "m01s30i462":"Vertical_integral_of_eastward_water_vapour_flux",
  "m01s30i463":"Vertical_integral_of_northward_water_vapour_flux",
  "STABLE BL INDICATOR":"Stable boundary layer indicator",
  "X-COMP SURFACE BL STRESS":"surface_downward_eastward_stress",
  "STRATOCUM. OVER STABLE BL INDICATOR":"Stratocumulus over stable boundary layer indicator",
  "WELL_MIXED BL INDICATOR":"Well-mixed boundary layer indicator",
  "Y-COMP SURFACE BL STRESS":"surface_downward_northward_stress",
  "DECOUPLED SC. OVER CU. INDICATOR":"Decoupled stratocumulus over cumulus indicator",
  "TOTAL CLOUD AMOUNT - RANDOM OVERLAP":"cloud_area_fraction_assuming_random_overlap",
  "TOTAL CLOUD AMOUNT MAX/RANDOM OVERLP":"cloud_area_fraction_assuming_maximum_random_overlap",
  "RADAR REFLECTIVITY AT 1KM AGL (dBZ)":"radar_reflectivity_due_to_all_hydrometeors_at_1km_altitude",
  "TURBULENT MIXING HT AFTER B.LAYER m":"Turbulent mixing height after boundary layer",
  "EVAP FROM OPEN SEA: SEA MEAN KG/M2/S":"Evaporation flux from open sea",
  "DECOUPLED SC. NOT OVER CU. INDICATOR":"Decoupled stratocumulus not over cumulus indicator",
  "CUMULUS-CAPPED BL INDICATOR":"Cumulus capped boundary layer indicator",
  }

units = {
  "m01s05i250":"Pa/s",
  "m01s05i251":"Pa/s",
  "m01s01i202":"W/m2",
  "m01s03i319":"m",
  "m01s03i359":"m",
  "m01s03i462":"m/s",
  "m01s03i465":"m/s",
  "m01s03i476":"1",
  "m01s04i113":"dBZ",
  "m01s04i114":"dBZ",
  "m01s04i115":"dBZ",
  "m01s04i116":"dBZ",
  "m01s04i117":"dBZ",
  "m01s09i202":"1",
  "m01s30i403":"kg/m2",
  "m01s30i461":"kg/m2",
  "m01s30i462":"kg/m/s",
  "m01s30i463":"kg/m/s",
  "m01s01i161":"K/hr",
  "m01s02i161":"K/hr",
}


def postprocess_output(date1,outname,stream):
    # postprocess data from crun starting on date1 into file outname.nc
    if stream == "pf" and outname=="radar_reflectivity":
      print("skipping radar reflectivity: not outputted for Mirai transect")
      return 
    cell_methods,variables = file_variables[outname] 
    # loop through times for file output
    if outname in ["cloud_fraction","density","mixing_ratios","specific_humidity","pressure_rho_grid","pressure_theta_grid","radar_reflectivity","potential_temperature","winds"]:
      orog=iris.load("/gws/nopw/j04/terramaris/emmah/coupled_2km/orog.pp")[0]
    for i in range(0,6,reinit_step_nc)[:1]:
      date2 = date1 + dt.timedelta(i)
      print(date2)
      times = [date2+ dt.timedelta(i/24) for i in np.arange(0,24*reinit_step_nc,reinit_step_pp)]
      files = [path+"/tm2a_%s%04d%02d%02d_%02d.nc"%(stream,t.year,t.month,t.day,t.hour)\
             for t in times]
      print(files)
#      if stream in ["pa","pc"]:
#         files.append(orogpath)
      data = iris.load(files[:6])
      for name in units.keys():
        for cube in data.extract(name):
          cube.units = units[name]
      for name in names.keys():
        for cube in data.extract(name):
            cube.rename(names[name])
      data = data.extract(variables).extract(cell_methods)
      data=data.concatenate()
      if outname in ["winds","specific_humidity"]:
        data = iris.cube.CubeList([cube for cube in data if cube.shape[1]==70])
      if outname == "surf_inst":
        data = iris.cube.CubeList([cube for cube in data if cube.shape[1]!=70])
      print(data)
      assert len(data) == len(variables),variables
      for cube in data:
        cube.coord("time").convert_units("hours since 2015-11-01 00:00:00")
        cube.var_name=None
        cube.remove_coord("latitude")
        cube.remove_coord("longitude")
        cube.coord("grid_latitude").rename("latitude")
        cube.coord("grid_longitude").rename("longitude")
        cube.coord("latitude").bounds  = None
        cube.coord("longitude").bounds  = None
        cube.coord("latitude").coord_system = iris.coord_systems.GeogCS(6371229.0)
        cube.coord("longitude").coord_system = iris.coord_systems.GeogCS(6371229.0)
        for c in cube.coords():
          c.var_name=None
        if cube.shape[1] == 70:
            orog2 = orog.regrid(cube,iris.analysis.Linear())
            orog_c2 = iris.coords.AuxCoord(orog2.data,long_name='surface_altitude',units='m')
            cube.add_aux_coord(orog_c2,[2,3])
            cube.coord("height_above_reference_ellipsoid").rename("level_height")
            cube.coord("Fraction of orographic height").rename("sigma")
            a = HybridHeightFactory(cube.coord('level_height'),\
                          cube.coord('sigma'),\
                          cube.coord('surface_altitude'))
            cube.add_aux_factory(a)
#            z = a.make_coord(cube.coord_dims)
#            cube.add_aux_coord(z,[1,2,3])
            #cube.coord("model_level_number").bounds=None
            #iris.util.promote_aux_coord_to_dim_coord(cube,"model_level_number")
            for name in ["level_height","sigma"]:#,"model_level_number"]:
              cube.coord(name).points = np.array(cube.coord(name).points)
   #   cube.coord(name).long_name = None
      iris.save(data,outpath+"%s/tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"
         %(stream,region[stream],outname,date2.year,date2.month,date2.day),zlib=True)


job =  int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
#for job in range(16,20):
var = list(file_variables.keys())[job]
print(var)
postprocess_output(dt.datetime(2015,12,1),var,"pf")
