#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH --array=[1-14]
#SBATCH -o /home/users/emmah/log/postproc_2km/pp1_%a.o
#SBATCH -e /home/users/emmah/log/postproc_2km/pp1_%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=24000

import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call

year=2015
t1 = dt.datetime(2015,12,1)

# Constraints on cell methods (ie UM time profiles)
c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="mean"))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not cube.cell_methods)

# system paths
path = "/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/%04d/%04d%02d%02dT0000Z/"%(year,t1.year,t1.month,t1.day)
orogpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/orog.pp"
outpath = "/work/scratch-nopw/emmah/postproc_2km_fc/"
finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cf309/"
# re-initialisation timesteps
reinit_step_pp = {"pa":4,"pb":24,"pc":4,"pd":12} # input, in hours
reinit_step_nc = {"pa":1,"pb":6,"pc":1,"pd":6}   # output, in days

# files and variables in the form:
# file name: (stream, [valid cell methods], [variable names])
file_variables = {
         "radar_reflectivity"   :("pa",[c_inst],["radar_reflectivity_due_to_all_hydrometeor_species"]),                        # 1
         "density"              :("pa",[c_inst],["air_density"]),                                                              # 2
         "potential_temperature":("pa",[c_inst],["air_potential_temperature"]),                                                # 3
         "pressure_rho_grid"    :("pa",[c_inst],["m01s00i407"]),                                                               # 4
         "pressure_theta_grid"  :("pa",[c_inst],["m01s00i408"]),                                                               # 5
         "cloud_fraction"       :("pa",[c_inst],["cloud_area_fraction_in_atmosphere_layer"]),                                  # 6
         "qcf"                  :("pa",[c_inst],["mass_fraction_of_cloud_ice_in_air"]),                                        # 7
         "qcl"                  :("pa",[c_inst],["mass_fraction_of_cloud_liquid_water_in_air"]),                               # 8
         "rain"                 :("pa",[c_inst],["mass_fraction_of_rain_in_air"]),                                             # 9
         "graupel"              :("pa",[c_inst],["mass_fraction_of_graupel_in_air"]),                                          #10
         "specific_humidity"    :("pa",[c_inst],["specific_humidity"]),                                                        #11
         "upward_air_velocity"  :("pa",[c_inst],["upward_air_velocity"]),                                                      #12
         "eastward_wind"        :("pa",[c_inst],["eastward_wind"]),                                                            #13
         "northward_wind"       :("pa",[c_inst],["northward_wind"]),                                                           #14
#
         "atmos":                ("pb",[c_mean,c_inst],["atmosphere_boundary_layer_thickness",                                 #15
                                        "toa_incoming_shortwave_flux",
                                        "toa_outgoing_shortwave_flux",
                                        "toa_outgoing_shortwave_flux_assuming_clear_sky",
                                        "toa_outgoing_longwave_flux",
                                        "toa_outgoing_longwave_flux_assuming_clear_sky",
                                        "Turbulent mixing height after boundary layer",
                                        "height_of_diagnostic_parcel_top",
                                        "combined_boundary_layer_type",
                                        "radar_reflectivity_due_to_all_hydrometeors_at_1km_altitude"]),
         "integral_inst":        ("pb",[c_inst],["Stable boundary layer indicator",                                            #16
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
         "integral_agg":         ("pb",[c_mean],["cloud_area_fraction",                                                        #17
                                        "very_low_type_cloud_area_fraction",
                                        "low_type_cloud_area_fraction",
                                        "medium_type_cloud_area_fraction",
                                        "high_type_cloud_area_fraction",
                                        "total_column_dry_mass",
                                        "atmosphere_mass_per_unit_area",
                                        "atmosphere_cloud_liquid_water_content",
                                        "atmosphere_cloud_ice_content",
                                        "Vertical_integral_of_eastward_water_vapour_flux",
                                        "Vertical_integral_of_northward_water_vapour_flux"]),
         "rainfall":             ("pb",[c_mean,c_inst],["stratiform_rainfall_amount",                                          #18
                                        "stratiform_snowfall_amount",
                                        "stratiform_rainfall_flux",
                                        "stratiform_snowfall_flux"]),
         "surf_inst":            ("pb",[c_inst],["surface_temperature",                                                        #19
                                        "surface_air_pressure",
                                        "eastward_wind",
                                        "northward_wind",
                                        "air_temperature",
                                        "specific_humidity",
                                        "relative_humidity",
                                        "dew_point_temperature",
                                        "canopy_height",
                                        "explicit_friction_velocity",
                                        "soil_moisture_content",
                                        "canopy_water_amount",
                                        "moisture_content_of_soil_layer",
                                        "air_pressure_at_sea_level"]),
         "surf_agg":             ("pb",[c_mean],["surface_net_downward_shortwave_flux",                                        #20
                                        "surface_net_downward_shortwave_flux",
                                        "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky",
                                        "surface_upwelling_shortwave_flux_in_air_assuming_clear_sky",
                                        "surface_downwelling_shortwave_flux_in_air",
                                        "surface_net_downward_longwave_flux",
                                        "surface_downwelling_longwave_flux_in_air",
                                        "surface_downwelling_longwave_flux_in_air_assuming_clear_sky",
                                        "surface_upward_sensible_heat_flux",
                                        "surface_upward_water_flux",
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
                                        "air_pressure_at_sea_level",
                                        "wind_speed_of_gust",
                                        "number_of_lightning_flashes"]),
         "range":             ("pb",[c_min,c_max],["air_temperature","hail_diameter"]),                                        #21
#                                         
         "dT_shortwave"          :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_shortwave_radiation_noPC2"]),    #22
         "dT_longwave"           :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_longwave_radiation_noPC2"]),     #23
         "dT_advection"          :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_advection"]),                    #24
         "dq_advection"          :("pc",[c_mean],["change_over_time_in_specific_humidity_due_to_advection"]),                  #25
         "w_T_rho"               :("pc",[c_mean],["w_T_rho_dry_1e-5"]),                                                        #26
         "w_q_rho"               :("pc",[c_mean],["w_q_rho_dry_1e-5"]),                                                        #27
         "dT_total"              :("pc",[c_mean],["change_over_time_in_air_temperature"]),                                     #28
         "dq_total"              :("pc",[c_mean],["change_over_time_in_specific_humidity"]),                                   #29
         "temperature_pc"        :("pc",[c_mean],["air_temperature"]),                                                         #30
         "specific_humidity_pc"  :("pc",[c_mean],["specific_humidity"]),                                                       #31
         "upward_air_velocity_pc":("pc",[c_mean],["m01s30i003"]),                                                              #32 
#
         "Heavyside"             :("pd",[c_inst],["Heaviside function on pressure levels"]),                                   #33
         "temperature_pd"        :("pd",[c_inst],["air_temperature"]),                                                         #34
         "geopotential_height"   :("pd",[c_inst],["geopotential_height"]),                                                     #35
         "omega"                 :("pd",[c_inst],["lagrangian_tendency_of_air_pressure"]),                                     #36
         "relative_humidity"     :("pd",[c_inst],["relative_humidity"]),                                                       #37
         "specific_humidity_pd"  :("pd",[c_inst],["specific_humidity"]),                                                       #38
         "theta_w"               :("pd",[c_inst],["wet_bulb_potential_temperature"]),                                          #39
         "eastward_wind_pd"      :("pd",[c_inst],["eastward_wind"]),                                                           #40
         "div_vort"              :("pd",[c_inst,c_mean],["divergence_of_wind","atmosphere_relative_vorticity"]),               #41
         "northward_wind_pd"     :("pd",[c_inst],["northward_wind"]),                                                          #42
         "potential_vorticity"   :("pd",[c_inst],["potential_vorticity"]),                                                     #43
         "daily_mean"            :("pd",[c_mean],["eastward_wind",                                                             #44
                                                  "northward_wind",
                                                  "wet_bulb_potential_temperature",
                                                  "Heaviside function on pressure levels"])
         
         }     


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
  "m01s30i034":"w_T_rho_dry_1e-5",
  "m01s30i035":"w_q_rho_dry_1e-5",
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
  "m01s30i034":"kg.K/m/s",
  "m01s30i035":"kg/m/s",
  "change_over_time_in_air_temperature_due_to_advection":"K/hr",
  "change_over_time_in_specific_humidity_due_to_advection":"1/hr",
  "change_over_time_in_air_temperature":"K/hr",
  "change_over_time_in_specific_humidity":"1/hr",
# pd
  "m01s15i229":"m2.K/s/kg",
  "divergence_of_wind":"10e6 s-1",
  "atmosphere_relative_vorticity":"10e6 s-1"
}


def postprocess_output(date1,outname):
    # postprocess data from crun starting on date1 into file outname.nc
    stream,cell_methods,variables = file_variables[outname] 
    # loop through times for file output
    for i in range(0,6,reinit_step_nc[stream])[:1]:
      date2 = date1 + dt.timedelta(i)
      print(date2)
      # times for file input
      times = [date2+ dt.timedelta(i/24) for i in np.arange(0,24*reinit_step_nc[stream],reinit_step_pp[stream])]
      files = [path+"/tm2a_%s%04d%02d%02d_%02d.pp"%(stream,t.year,t.month,t.day,t.hour)\
             for t in times]
      if stream in ["pa","pc"]:
         # ensures 3d altitude aux coordinate will be build
         files.append(orogpath)
      data = iris.load(files[0])
      for cube in data:
        if cube.name() in units.keys():
            # correct units
            cube.units = units[cube.name()]
        if cube.name() in names.keys():
            # correct variable names
            cube.rename(names[cube.name()])
      data=data.concatenate()
      data=data.merge()
      if stream=="pd" and cube.name()!="Heaviside function on pressure levels":
         # extract heaviside function for masking data below orography
         H = data.extract("Heaviside function on pressure levels")[0].copy()
      # extract variables for saving 
      data = data.extract(variables).extract(cell_methods)
      print(data)
      # check for the right number of variables. "range" has 2x "variables" due to max and mins
      if outname == "range": 
        assert len(data) == len(variables)*2
      else:
        assert len(data) == len(variables)
      # mask pressure level data below orography
      if outname =="daily_mean":
        if len(data.extract("Heaviside function on pressure levels"))>0: 
          H = data.extract("Heaviside function on pressure levels")[0]
          for cube in data:
            if cube.name() not in ["Heaviside function on pressure levels","wet_bulb_potential_temperature"]: # wbpt is on a different grid
              cube.data.mask += (H.data==0)
              #cube.data = cube.data/H.data
        # if heaviside wasnt averaged, approximate it
        else:
          from iris.coord_categorisation import add_day_of_year
          add_day_of_year(H,"time","doyr")
          H=H.aggregated_by("doyr",iris.analysis.MEAN)[:-1]
          for cube in data:
            if cube.name() != "wet_bulb_potential_temperature": # wbpt is on a different grid
              cube.data = np.ma.masked_array(cube.data,mask=(H.data==0))
#              cube.data = cube.data/H.data
      elif outname == "div_vort":
        cp = iris.Constraint(pressure=lambda p: p in  [200,700,850,925])   # div and vort are available on a smaller number of levels
        H = H.extract(cp)
        data[0].data = np.ma.masked_array(data[0].data,mask=1-H.data)
        data[1].data = np.ma.masked_array(data[1].data,mask=1-H.data)
      elif stream=="pd" and outname != "Heaviside" and outname != "theta_w":
        data[0].data = np.ma.masked_array(data[0].data,mask=1-H.data)
      for cube in data:
        # convert units to hour since simulation started
        cube.coord("time").convert_units("hours since 2015-11-01 00:00:00")
      # save to netcdf
      outfile = "tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"%(stream,outname,date2.year,date2.month,date2.day)
      iris.save(data,outpath+"%s/%s"%(stream,outfile),zlib=True)

def copy_remove(date1,outname):
    # postprocess data from crun starting on date1 into file outname.nc
    stream,cell_methods,variables = file_variables[outname] 
    # loop through times for file output
    for i in range(0,6,reinit_step_nc[stream]):
      date2 = date1 + dt.timedelta(i)
      outfile = "tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"%(stream,outname,date2.year,date2.month,date2.day)
      cmd = "cp %s/%s/%s %s/%s/%s"%(outpath,stream,outfile,finalpath,stream,outfile)
      check_call(cmd,shell=True)
      cmd = "rm %s/%s/%s"%(outpath,stream,outfile)
      check_call(cmd,shell=True)


job =  int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
#for job in range(31,32):
#for job in range(16,20):
var = list(file_variables.keys())[job]
print(var)
#if var in ["Heavyside","temperature_pd","geopotential_height","omega","relative_humidity","theta_w","updraft_mass_flux"]:
#  exit()
postprocess_output(t1,var)
#copy_remove(t1,var)










