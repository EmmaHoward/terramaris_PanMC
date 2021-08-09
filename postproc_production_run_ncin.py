#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH --array=[21]
#SBATCH -o /home/users/emmah/log/postproc_N1280/pp2_%a.o
#SBATCH -e /home/users/emmah/log/postproc_N1280/pp2_%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=24000

import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call
from iris.aux_factory import HybridHeightFactory
year=2015
t1 = dt.datetime(2015,12,7)

# Constraints on cell methods (ie UM time profiles)
c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method  in ["mean","sum"]))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not (cube.cell_methods) or (cube.cell_methods[0].method=="point"))

# system paths
path = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/u-cf309/%04d/%04d%02d%02dT0000Z/"%(year,t1.year,t1.month,t1.day)
orogpath = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"
outpath = "/work/scratch-nopw/emmah/postproc_N1280_fc/"
finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/production_runs/201516_u-cf309/"
# re-initialisation timesteps
reinit_step_pp = {"pa":4,"pb":24,"pc":6,"pd":12} # input, in hours
reinit_step_nc = {"pa":1,"pb":6,"pc":1,"pd":6}   # output, in days

# files and variables in the form:
# file name: (stream, [valid cell methods], [variable names])
file_variables = {
         "radar_reflectivity"   :("pa",[c_inst],["radar_reflectivity_due_to_all_hydrometeor_species"]),                        # 1
         "updraft_mass_flux"    :("pa",[c_mean],["atmosphere_updraft_convective_mass_flux"]),                                  # 2
         "downdraft_mass_flux"  :("pa",[c_mean],["atmosphere_downdraft_convective_mass_flux"]),                                # 3
         "density"              :("pa",[c_inst],["air_density"]),                                                              # 4
         "potential_temperature":("pa",[c_inst],["air_potential_temperature"]),                                                # 5
         "pressure_rho_grid"    :("pa",[c_inst],["m01s00i407"]),                                                               # 6
         "pressure_theta_grid"  :("pa",[c_inst],["m01s00i408"]),                                                               # 7
         "cloud_fraction"       :("pa",[c_inst],["cloud_area_fraction_in_atmosphere_layer"]),                                  # 8
         "qcf"                  :("pa",[c_inst],["mass_fraction_of_cloud_ice_in_air"]),                                        # 9
         "qcl"                  :("pa",[c_inst],["mass_fraction_of_cloud_liquid_water_in_air"]),                               #10
         "rain"                 :("pa",[c_inst],["mass_fraction_of_rain_in_air"]),                                             #11
         "specific_humidity"    :("pa",[c_inst],["specific_humidity"]),                                                        #12
         "upward_air_velocity"  :("pa",[c_inst],["upward_air_velocity"]),                                                      #13
         "eastward_wind"        :("pa",[c_inst],["eastward_wind"]),                                                            #14
         "northward_wind"       :("pa",[c_inst],["northward_wind"]),                                                           #15
#
         "atmos":                ("pb",[c_mean,c_inst],["atmosphere_boundary_layer_thickness",                                 #16
                                        "toa_incoming_shortwave_flux",
                                        "toa_outgoing_shortwave_flux",
                                        "toa_outgoing_shortwave_flux_assuming_clear_sky",
                                        "toa_outgoing_longwave_flux",
                                        "toa_outgoing_longwave_flux_assuming_clear_sky",
                                        "Turbulent mixing height after boundary layer",
                                        "height_of_diagnostic_parcel_top",
                                        "combined_boundary_layer_type",
                                        "radar_reflectivity_due_to_all_hydrometeors_at_1km_altitude"]),
         "integral_inst":        ("pb",[c_inst],["Stable boundary layer indicator",                                            #17
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
         "integral_agg":         ("pb",[c_mean,c_inst],["cloud_area_fraction",                                                        #18
                                        "very_low_type_cloud_area_fraction",
                                        "low_type_cloud_area_fraction",
                                        "medium_type_cloud_area_fraction",
                                        "high_type_cloud_area_fraction",
                                        #"total_column_dry_mass",
                                        #"atmosphere_mass_per_unit_area",
                                        #"atmosphere_cloud_liquid_water_content",
                                        #"atmosphere_cloud_ice_content",
                                        "Vertical_integral_of_eastward_water_vapour_flux",
                                        "Vertical_integral_of_northward_water_vapour_flux"]),
         "rainfall":                 ("pb",[c_mean,c_inst],["stratiform_rainfall_amount",                                     #19
                                        "stratiform_snowfall_amount",
                                        "stratiform_rainfall_flux",
                                        "stratiform_snowfall_flux",
                                        "convective_rainfall_amount",
                                        "convective_snowfall_amount",
                                        "convective_rainfall_flux",
                                        "convective_snowfall_flux"]),
         "surf_inst":            ("pb",[c_inst],["surface_temperature",                                                        #20
                                        "surface_air_pressure",
                                        "eastward_wind",
                                        "northward_wind",
                                        "air_temperature",
                                        "specific_humidity",
                                        "relative_humidity",
                                        "stomatal_conductance", 
                                        "dew_point_temperature",
                                        "canopy_height",
                                        "explicit_friction_velocity",
                                        "soil_moisture_content",
                                        "canopy_water_amount",
                                        "moisture_content_of_soil_layer",
                                        "air_pressure_at_sea_level"]),
         "surf_agg":             ("pb",[c_mean],["m01s01i201",                                                                 #21
                                        "m01s01i202",
                                        "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky",
                                        "surface_upwelling_shortwave_flux_in_air_assuming_clear_sky",
                                        "surface_downwelling_shortwave_flux_in_air",
                                        "surface_net_downward_longwave_flux",
                                        "surface_downwelling_longwave_flux_in_air",
                                        "surface_downwelling_longwave_flux_in_air_assuming_clear_sky",
                                        "surface_upward_sensible_heat_flux",
                                        "surface_upward_water_flux",
                                        "surface_moisture_flux", # 1 only
                                        "eastward_wind",
                                        "northward_wind",
                                        "wind_speed",
                                        "Evaporation flux from open sea",
                                        "surface_upward_latent_heat_flux",
                                        "air_temperature",
                                        "specific_humidity",
                                        "relative_humidity",
                                        "Evaporation from soil surface",
                                        "Evaporation from canopy",
                                        "surface_downward_eastward_stress",
                                        "surface_downward_northward_stress",
                                        "transpiration_amount",
                                        "surface_runoff_flux",
                                        "subsurface_runoff_flux",
                                        "air_pressure_at_sea_level"]),
         "range":             ("pb",[c_min,c_max],["air_temperature","wind_speed_of_gust"]),                                                         #22

         "dT_shortwave"	         :("pc",[c_inst],["m01s01i161"]),                                                               #23
         "dT_longwave"           :("pc",[c_inst],["m01s02i161"]),                                                               #24
         "dT_advection"          :("pc",[c_inst],["tendency_of_air_temperature_due_to_advection"]),                             #25
#         "dT_advection"          :("pc",["change_over_time_in_air_temperature_due_to_advection"]),             
         "dq_advection"          :("pc",[c_inst],["tendency_of_specific_humidity_due_to_advection"]),                           #26
#         "dq_advection"          :("pc",["change_over_time_in_specific_humidity_due_to_advection"]), 
         "w_T_rho"               :("pc",[c_inst],["m01s30i034"]),                                                               #27
         "w_q_rho"               :("pc",[c_inst],["m01s30i035"]),                                                               #28
         "dT_total"              :("pc",[c_inst],["tendency_of_air_temperature"]),                                              #29
#         "dT_total"              :("pc",["change_over_time_in_air_temperature"]),                      
         "dq_total"              :("pc",[c_inst],["tendency_of_specific_humidity"]),                                            #30
#         "dq_total"              :("pc",["change_over_time_in_specific_humidity"]),                    
         "temperature_pc"        :("pc",[c_inst],["air_temperature"]),                                                          #31
         "specific_humidity_pc"  :("pc",[c_inst],["specific_humidity"]),                                                        #32
         "upward_air_velocity_pc":("pc",[c_inst],["m01s30i003"]),#"upward_air_velocity"]),                                      #33
#         
         "Heaviside"             :("pd",[c_inst],["HEAVYSIDE FN ON P LEV/UV GRID"]),                                        #34
         "temperature_pd"        :("pd",[c_inst],["air_temperature"]),                                                      #35
         "geopotential_height"   :("pd",[c_inst],["geopotential_height"]),                                                  #36
         "omega"                 :("pd",[c_inst],["lagrangian_tendency_of_air_pressure"]),                                  #37
         "relative_humidity"     :("pd",[c_inst],["relative_humidity"]),                                                    #38
         "specific_humidity_pd"  :("pd",[c_inst],["specific_humidity"]),                                                    #39
         "theta_w"               :("pd",[c_inst],["wet_bulb_potential_temperature"]),                                       #40
         "eastward_wind_pd"      :("pd",[c_inst],["eastward_wind"]),                                                        #41
         "northward_wind_pd"     :("pd",[c_inst],["northward_wind"])}                                                       #42




names = {
# multiple
  "x_wind":"eastward_wind",
  "y_wind":"northward_wind",
# pa
  "TOTAL RADAR REFLECTIVITY 3D":"radar_reflectivity_due_to_all_hydrometeor_species",
#  "m01s05i250":"atmosphere_updraft_convective_mass_flux",
#  "m01s05i251":"atmosphere_downdraft_convective_mass_flux",
# pb
  "m01s01i202":"surface_net_downward_shortwave_flux",
  "m01s03i241":"surface_upward_water_flux",
  "m01s03i319":"canopy_height",
  "m01s03i359":"height_of_diagnostic_parcel_top",
  "m01s03i462":"stomatal_conductance",
  "m01s03i465":"explicit_friction_velocity",
  "m01s03i476":"combined_boundary_layer_type",
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
# pc
  "m01s01i161":"change_over_time_in_air_temperature_due_to_shortwave_radiation_noPC2",
  "m01s02i161":"change_over_time_in_air_temperature_due_to_longwave_radiation_noPC2",
  "m01s30i034":"w_T_rho_dry_1e-5",
  "m01s30i035":"w_q_rho_dry_1e-5",
# pd
  "HEAVYSIDE FN ON P LEV/UV GRID":"Heavyside function on pressure levels"

  }

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
  "m01s09i202":"1",
  "m01s30i403":"kg/m2",
  "m01s30i461":"kg/m2",
  "m01s30i462":"kg/m/s",
  "m01s30i463":"kg/m/s",
# pc
  "m01s01i161":"K/hr",
  "m01s02i161":"K/hr",
  "m01s30i034":"kg.K/m/s",
  "m01s30i035":"kg/m/s",
  "change_over_time_in_air_temperature_due_to_advection":"K/hr",
  "change_over_time_in_specific_humidity_due_to_advection":"1/hr",
  "change_over_time_in_air_temperature":"K/hr",
  "change_over_time_in_specific_humidity":"1/hr"
}


def postprocess_output(date1,outname):
    # postprocess data from crun starting on date1 into file outname.nc
    stream,cell_methods,variables = file_variables[outname] 
    # loop through times for file output
    if stream in ["pa","pc"]:
      orog=iris.load("/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp")[0]
      orog_c = iris.coords.AuxCoord(orog.data,long_name='surface_altitude',units='m')
    for i in range(0,6,reinit_step_nc[stream]):
      date2 = date1 + dt.timedelta(i)
      print(date2)
      times = [date2+ dt.timedelta(i/24) for i in np.arange(0,24*reinit_step_nc[stream],reinit_step_pp[stream])]
      files = [path+"/tma_%s%04d%02d%02d_%02d.nc"%(stream,t.year,t.month,t.day,t.hour)\
             for t in times]
#      if stream in ["pa","pc"]:
#         files.append(orogpath)
      data = iris.load(files)
      for name in units.keys():
        for cube in data.extract(name):
          cube.units = units[name]
      for name in names.keys():
        for cube in data.extract(name):
            cube.rename(names[name])
      data = data.extract(variables).extract(cell_methods)
      data=data.concatenate()
      print(data)
      #assert len(data) == len(variables),variables
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
        if stream=="pd" and cube.name()!="Heaviside function on pressure levels":
           heaviside = iris.load(files,"HEAVYSIDE FN ON P LEV/UV GRID")
           heaviside = heaviside.concatenate_cube()
           cube.data.mask += 1-heaviside.data
        for c in cube.coords():
          c.var_name=None
        if stream in ["pa","pc"]:
            if cube.name() in ["eastward_wind","northward_wind"]:
              orog2 = orog.regrid(cube,iris.analysis.Linear())
              orog_c2 = iris.coords.AuxCoord(orog2.data,long_name='surface_altitude',units='m')
              cube.add_aux_coord(orog_c2,[2,3])
            else:
              cube.add_aux_coord(orog_c,[2,3])
            a = HybridHeightFactory(cube.coord('height_above_reference_ellipsoid'),\
                          cube.coord('Fraction of orographic height'),\
                          cube.coord('surface_altitude'))
            z = a.make_coord(cube.coord_dims)
#            cube.add_aux_coord(z,[1,2,3])
            cube.coord("height_above_reference_ellipsoid").rename("level_height")
            cube.coord("Fraction of orographic height").rename("sigma")
            cube.add_aux_coord(z,[1,2,3])
            #cube.coord("model_level_number").bounds=None
            #iris.util.promote_aux_coord_to_dim_coord(cube,"model_level_number")
            for name in ["level_height","sigma"]:#,"model_level_number"]:
              cube.coord(name).points = np.array(cube.coord(name).points)
   #   cube.coord(name).long_name = None
      iris.save(data,outpath+"%s/tma_N1280_KPPcoupled_%s_%s_%04d%02d%02d.nc"
         %(stream,stream,outname,date2.year,date2.month,date2.day),zlib=True)


job =  int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
#for job in range(16,20):
var = list(file_variables.keys())[job]
print(var)
postprocess_output(dt.datetime(2015,12,7),var)
