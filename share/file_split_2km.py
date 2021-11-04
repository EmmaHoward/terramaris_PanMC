
import iris

# Constraints on cell methods (ie UM time profiles)
c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method  in ["mean","sum"]))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not (cube.cell_methods) or (cube.cell_methods[0].method=="point"))

packing = {
"00002":-8,
"00003":-8,
"00004":-8,
"00010":-24,
"00150":-10,
"00389":-20,
"00407":-3,
"00408":-3,
"00265":-6,
"00266":-6,
"00267":-6,
"00268":-6,
"04100":-99,
"00012":-24,
"00254":-24,
"00272":-24,
"00273":-24,
"04113":-3,
"04114":-3,
"04115":-3,
"04116":-3,
"04117":-3,
"00025":-3,
"01207":-6,
"01208":-6,
"01209":-6,
"02205":-6,
"02206":-6,
"03304":-3,
"03359":-3,
"03476":-10,
"09202":-6,
"09203":-6,
"09204":-6,
"09205":-6,
"30462":-3,
"30463":-3,
"03305":-3,
"03306":-3,
"03307":-3,
"03308":-3,
"03309":-3,
"03310":-3,
"09216":-6,
"09217":-6,
"30403":-3,
"30404":-3,
"30405":-8,
"30406":-8,
"30461":-8,
"04201":-8,
"04202":-8,
"04203":-18,
"04204":-18,
"00024":-8,
"00409":-3,
"03225":-6,
"03226":-6,
"03236":-6,
"03237":-24,
"03245":-6,
"03250":-6,
"03465":-10,
"01201":-12,
"01202":-6,
"01210":-99,
"01211":-99,
"01235":-12,
"02201":-12,
"02207":-12,
"02208":-99,
"03217":-12,
"03232":-24,
"03234":-12,
"03241":-12,
"03460":-10,
"03461":-10,
"03463":-99,
"21104":-1,
   }


file_variables_pepf = {
         "cloud_fraction"       :([c_inst],["cloud_area_fraction_in_atmosphere_layer",                               # 1
                                         "cloud_volume_fraction_in_atmosphere_layer",
                                         "liquid_cloud_volume_fraction_in_atmosphere_layer",
                                         "ice_cloud_volume_fraction_in_atmosphere_layer"]),
#                                         "ice_aggregate_fraction"]),

         "radar_reflectivity"   :([c_inst],["radar_reflectivity_due_to_graupel_alone",                               # 2
                                          "radar_reflectivity_due_to_ice_alone",
#                                          "radar_reflectivity_due_to_cloud_snow_alone",
                                          "radar_reflectivity_due_to_rain_alone",
                                          "radar_reflectivity_due_to_cloud_liquid_alone"]),
         "mixing_ratios"        :([c_inst],["mass_fraction_of_cloud_ice_in_air",                                     # 3
                                                 "mass_fraction_of_cloud_liquid_water_in_air",
                                                 "mass_fraction_of_graupel_in_air",
                                                 "mass_fraction_of_rain_in_air"]),
         "density"              :([c_inst],["air_density"]),                                                         # 4
         "potential_temperature":([c_inst],["air_potential_temperature"]),                                           # 5
         "pressure_rho_grid"    :([c_inst],["m01s00i407"]),                                                          # 6
         "pressure_theta_grid"  :([c_inst],["m01s00i408"]),                                                          # 7
         "specific_humidity"    :([c_inst],["specific_humidity"]),                                                   # 8
         "eastward_wind"        :([c_inst],["eastward_wind"]),                                                       # 9
         "northward_wind"       :([c_inst],["northward_wind"]),                                                      #10
         "upward_air_velocity"  :([c_inst],["upward_air_velocity"]),                                                 #11
#
         "atmos":                ([c_mean,c_inst],["atmosphere_boundary_layer_thickness",                            #12
                                        "toa_incoming_shortwave_flux",
                                        "toa_outgoing_shortwave_flux",
                                        "toa_outgoing_shortwave_flux_assuming_clear_sky",
                                        "toa_outgoing_longwave_flux",
                                        "toa_outgoing_longwave_flux_assuming_clear_sky",
                                        "Turbulent mixing height after boundary layer",
                                        "height_of_diagnostic_parcel_top",
                                        "combined_boundary_layer_type"]),
         "integral_inst":        ([c_inst],["Stable boundary layer indicator",                                       #13
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
         "integral_agg":         ([c_mean],[                                                                         #14
                                        "very_low_type_cloud_area_fraction",
                                        "low_type_cloud_area_fraction",
                                        "medium_type_cloud_area_fraction",
                                        "high_type_cloud_area_fraction",
                                        "Vertical_integral_of_eastward_water_vapour_flux",
                                        "Vertical_integral_of_northward_water_vapour_flux"]),
         "rainfall":                 ([c_mean,c_inst],["stratiform_rainfall_amount",                                 #15
                                        "stratiform_rainfall_flux",
                                        "stratiform_snowfall_amount",
                                        "stratiform_snowfall_flux"]),
         "surf_inst":            ([c_inst],["surface_temperature",                                                   #16
                                        "surface_air_pressure",
                                        "m01s03i225",
                                        "m01s03i226",
                                        "air_temperature",
                                        "specific_humidity",
                                        "relative_humidity",
                                        "dew_point_temperature",
                                        "explicit_friction_velocity"]),
         "surf_agg":             ([c_mean,c_max,c_inst],[
                                        "m01s01i201",                                                                #17
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

file_variables_full = {
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
                                        "m01s03i225",
                                        "m01s03i226",
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
         "surf_agg":             ("pb",[c_mean],["m01s01i201",                                                                 #20
                                        "m01s01i202",
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
                                        "Number_of_lightning_flashes"]),
         "range":             ("pb",[c_min,c_max],["air_temperature","maximum_predicted_hailstone_size_at_surface","wind_speed_of_gust"]),                   #21
#                                         
         "dT_shortwave"          :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_shortwave_radiation_noPC2"]),    #22
         "dT_longwave"           :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_longwave_radiation_noPC2"]),     #23
         "dT_advection"          :("pc",[c_mean],["change_over_time_in_air_temperature_due_to_advection"]),                    #24
         "dq_advection"          :("pc",[c_mean],["change_over_time_in_specific_humidity_due_to_advection"]),                  #25
         "w_T_rho"               :("pc",[c_mean],["w_T_rho_1e-5"]),                                                        #26
         "w_q_rho"               :("pc",[c_mean],["w_q_rho_1e-5"]),                                                        #27
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
#         "potential_vorticity"   :("pd",[c_inst],["potential_vorticity"]),                                                     #
         "daily_mean"            :("pd",[c_mean],["eastward_wind",                                                             #43
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


