#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python

import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call
import sys


# Constraints on cell methods (ie UM time profiles)
c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="mean"))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not cube.cell_methods)

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
                                        "number_of_turbulent_mixing_levels",
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
         "integral_agg":         ("pb",[c_mean],["cloud_area_fraction",                                                        #18
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
                                        "dew_point_temperature",
                                        "canopy_height",
                                        "explicit_friction_velocity",
                                        "soil_moisture_content",
                                        "canopy_water_amount",
                                        "stomatal_conductance", 
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
         "range":             ("pb",[c_min,c_max],["air_temperature","wind_speed_of_gust"]),                                   #22
#                                         
         "dT_shortwave"          :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_shortwave_radiation_noPC2"]),    #23
         "dT_longwave"           :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_longwave_radiation_noPC2"]),     #24
         "dT_advection"          :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_advection"]),                    #25
         "dq_advection"          :("pc",[c_mean],["change_over_time_in_specific_humidity_due_to_advection"]),                  #26
        "w_T_rho"               :("pc",[c_mean],["w_T_rho_dry_1e-5"]),                                                        #27
         "w_q_rho"               :("pc",[c_mean],["w_q_rho_dry_1e-5"]),                                                        #28
         "dT_total"              :("pc",[c_mean],["change_over_time_in_air_temperature"]),                                     #29
         "dq_total"              :("pc",[c_mean],["change_over_time_in_specific_humidity"]),                                   #30
         "temperature_pc"        :("pc",[c_mean],["air_temperature"]),                                                         #31
         "specific_humidity_pc"  :("pc",[c_mean],["specific_humidity"]),                                                       #32
         "upward_air_velocity_pc":("pc",[c_mean],["m01s30i003"]),                                                              #33 
#
         "Heavyside"             :("pd",[c_inst],["Heaviside function on pressure levels"]),                                   #34
         "temperature_pd"        :("pd",[c_inst],["air_temperature"]),                                                         #35
         "geopotential_height"   :("pd",[c_inst],["geopotential_height"]),                                                     #36
         "omega"                 :("pd",[c_inst],["lagrangian_tendency_of_air_pressure"]),                                     #37
         "relative_humidity"     :("pd",[c_inst],["relative_humidity"]),                                                       #38
         "specific_humidity_pd"  :("pd",[c_inst],["specific_humidity"]),                                                       #39
         "theta_w"               :("pd",[c_inst],["wet_bulb_potential_temperature"]),                                          #40
         "eastward_wind_pd"      :("pd",[c_inst],["eastward_wind"]),                                                           #41
         "div_vort"              :("pd",[c_inst,c_mean],["divergence_of_wind","atmosphere_relative_vorticity"]),               #42
         "northward_wind_pd"     :("pd",[c_inst],["northward_wind"]),                                                          #43
         "potential_vorticity"   :("pd",[c_inst],["potential_vorticity"]),                                                     #44
         "daily_mean"            :("pd",[c_mean],["eastward_wind",                                                             #45
                                                  "northward_wind",
                                                  "wet_bulb_potential_temperature",
                                                  "Heaviside function on pressure levels"])
         
         }     



file_variables_ocean = {
         "sea_water_temperature":["sea_water_temperature"],                          # 1
         "sea_water_salinity"   :["sea_water_salinity"],                             # 2
         "heat_flux_correction" :["heat_flux_correction_per_metre"],                 # 3
         "salinity_correction"  :["virtual_salt_flux_correction_per_metre"],         # 4
         "wT_turb"              :["Turbulent_vertical_sea_water_temperature_flux"],  # 5
         "wT_solar"             :["Radiative_vertical_sea_water_temperature_flux"],  # 6
         "wS_turb"              :["Turbulent_vertical_salinity_flux"],               # 7
         "sea_water_density"    :["sea_water_density"],                              # 8
         "cp"                   :["specific_heat_capacity_of_sea_water"],            # 9
         "mixedlayerdepth"      :["ocean_mixed_layer_thickness"],                    #10
         "sea_surface"          :["surface_downward_eastward_stress",
                                  "surface_downward_northward_stress",
                                  "surface_net_downward_shortwave_flux",
                                  "surface_downward_heat_flux_in_sea_water_excludes_shortwave",
                                  "precipitation minus evaporation"]}


