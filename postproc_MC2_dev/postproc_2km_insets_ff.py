#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
import sys
import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call
from iris.aux_factory import HybridHeightFactory

print(sys.argv[1])
i_start = int(sys.argv[2])
i_end = int(sys.argv[3])
stream = sys.argv[4]

t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1


region = {"pe":"Java","pf":"Bengkulu"}
#region = {"ne":"Java","nf":"Bengkulu"}

# Constraints on cell methods (ie UM time profiles)
c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method  in ["mean","sum"]))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not (cube.cell_methods) or (cube.cell_methods[0].method=="point"))

from file_split_2km import file_variables_pepf as file_variables
# system paths
orogpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/orog.pp"
outpath = "/work/scratch-nopw/emmah/postproc_2km_fc/"
finalpath = "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/"
# re-initialisation timesteps
#if t1.day in [1,2]:
#  reinit_step_pp = {"pa":4,"pb":4,"pc":4,"pd":6} # input, in hours
#elif t1.day in [3,4,5,6]:
#  reinit_step_pp = {"pa":1,"pb":4,"pc":1,"pd":6} # input, in hours
#reinit_step_nc = {"pa":1,"pb":6,"pc":1,"pd":6}   # output, in days
#reinit_step_nc = {"pa":4,"pb":24,"pc":4,"pd":24}   # output, in days

reinit_step_nc=1

names = {
# multiple
  "x_wind":"eastward_wind",
  "y_wind":"northward_wind",
  "m01s01i202":"surface_net_downward_shortwave_flux",
  "m01s03i241":"surface_upward_water_flux",
  "m01s03i359":"height_of_diagnostic_parcel_top",
  "m01s03i465":"explicit_friction_velocity",
  "m01s03i476":"combined_boundary_layer_type",
  "m01s09i202":"very_low_type_cloud_area_fraction",
  "m01s30i403":"total_column_dry_mass",
  "m01s30i461":"atmosphere_mass_content_of_water_vapor",
  "m01s30i462":"Vertical_integral_of_eastward_water_vapour_flux",
  "m01s30i463":"Vertical_integral_of_northward_water_vapour_flux",
  "Number_of_lightning_flashes":"number_of_lightning_flashes",
  "radar_reflectivity_due_to_ice_aggregates_alone":"radar_reflectivity_due_to_ice_alone",
  "radar_reflectivity_due_to_cloud_alone":"radar_reflectivity_due_to_cloud_liquid_alone",
  }


units = {
  "m01s01i202":"W/m2",
  "m01s03i359":"m",
  "m01s03i465":"m/s",
  "m01s03i476":"1",
  "m01s09i202":"1",
  "m01s30i403":"kg/m2",
  "m01s30i461":"kg/m2",
  "m01s30i462":"kg/m/s",
  "m01s30i463":"kg/m/s",
}


def postprocess_output(date1,outname,stream):
    path = "/work/scratch-nopw/emmah/pepf/%04d%02d%02dT0000Z/"%(date1.year,date1.month,date1.day)
    # postprocess data from crun starting on date1 into file outname.nc
    #if stream == "pf" and outname=="radar_reflectivity":
    #  print("skipping radar reflectivity: not outputted for Mirai transect")
    #  return 
    cell_methods,variables = file_variables[outname] 
    # loop through times for file output
    if i_start==0 and i_end==6:
      hourly=False
      index_range =range(0,6*24,24)
    else:
      hourly=True
      index_range = range(i_start*24,i_end*24,1) 
    for i in index_range:
      date2 = date1 + dt.timedelta(i/24)
#     if not os.path.exists(outpath+"%s/tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"
#         %(stream,region[stream],outname,date2.year,date2.month,date2.day)):
      print(date2)
#      times = [date2+ dt.timedelta(i/24/60) for i in np.arange(1,60*24*reinit_step_nc,reinit_step_pp)]
#      times[0] = times[0].replace(minute=0)
      if hourly:
        files = path+"/tm2a_%s%04d%02d%02d_%02d??"%(stream,date2.year,date2.month,date2.day,date2.hour)
      else:
        files = path+"/tm2a_%s%04d%02d%02d_????"%(stream,date2.year,date2.month,date2.day)
#      files = [path+"/tm2a_%s%04d%02d%02d_%02d%02d"%(stream,t.year,t.month,t.day,t.hour,t.minute)\
#             for t in times][:10]
#      print(files)
#      if stream in ["pa","pc"]:
#         files.append(orogpath)
      data2 = iris.load(files)#,callback=callback)
      for name in units.keys():
        for cube in data2.extract(name):
          cube.units = units[name]
      for name in names.keys():
        for cube in data2.extract(name):
            cube.rename(names[name])
      data = data2.extract(variables).extract(cell_methods)
      data=data.concatenate()
      if outname in ["winds","specific_humidity"]:
        data = iris.cube.CubeList([cube for cube in data if cube.shape[1]==70])
      if outname == "surf_inst":
        data = iris.cube.CubeList([cube for cube in data if cube.shape[1]!=70])
#      assert len(data) == len(variables),variables
      print(data)
      for cube in data:
        # convert units to hour since simulation started
        cube.coord("time").convert_units("hours since 2015-11-01 00:00:00")
        # convert to float # already done on Archer2
#        pack = packing["%02d%03d"%(cube.attributes["STASH"].section,cube.attributes["STASH"].item)]
#        if pack >-99:
#          cube.data = (np.round(cube.data/(2**pack))*2**pack).astype("float32")
#
#        if date1.day in [-1]:# netcdf files
#          cube.var_name=None
#          cube.long_name=None
#          cube.remove_coord("latitude")
#          cube.remove_coord("longitude")
#          cube.coord("grid_latitude").rename("latitude")
#          cube.coord("grid_longitude").rename("longitude")
#          cube.coord("latitude").bounds  = None
#          cube.coord("longitude").bounds  = None
#          cube.coord("latitude").coord_system = iris.coord_systems.GeogCS(6371229.0)
#          cube.coord("longitude").coord_system = iris.coord_systems.GeogCS(6371229.0)
#          for c in cube.coords():
#            c.var_name=None
#          if cube.shape[1] == 70:
#            cube.remove_coord("atmosphere_hybrid_height_coordinate")
#            cube.coord("model_level_number").long_name=None
#            cube.coord("height_above_reference_ellipsoid").rename("level_height")
#            cube.coord("Fraction of orographic height").rename("sigma")
#            cube.coord("model_level_number").bounds=None
#            iris.util.promote_aux_coord_to_dim_coord(cube,"model_level_number")
#            for name in ["level_height","sigma"]:#,"model_level_number"]:
#              cube.coord(name).points = np.array(cube.coord(name).points)
#        else:
        cube.remove_coord("forecast_period")
        cube.remove_coord("forecast_reference_time")
        if len(cube.cell_methods) > 0:
          if cube.cell_methods[0].method  in ["mean","sum"]:
            cube.cell_methods = (iris.coords.CellMethod("mean","time","5 minutes"),)
          elif cube.cell_methods[0].method == "point":
            cube.cell_methods = ()
          elif cube.cell_methods[0].method == "maximum":
            cube.cell_methods = (iris.coords.CellMethod("maximum","time","5 minutes"),)
        for c in cube.coords():
          if c.points.dtype in [int,"int32"]:
            c.points = c.points.astype("int32")
            if c.has_bounds():
              c.bounds = c.bounds.astype("int32")
          else:
            c.points = c.points.astype("float32")
            if c.has_bounds():
              c.bounds = c.bounds.astype("float32")         
          if c.name() in ["latitude","longitude"]:
            c.points = np.round(c.points * 100)/100
      if hourly:
        iris.save(data,outpath+"%s/tma_2km_KPPcoupled_%s_%s_%04d%02d%02d_%02d.nc"
         %(stream,region[stream],outname,date2.year,date2.month,date2.day,date2.hour),zlib=True)
      else:
        iris.save(data,outpath+"%s/tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"
         %(stream,region[stream],outname,date2.year,date2.month,date2.day),zlib=True)
 
job = int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
#for job in range(16,20):
var = list(file_variables.keys())[job]
print(var)

postprocess_output(t1,var,stream)
#postprocess_output(dt.datetime(2015,12,2),var,"pe")
#postprocess_output(dt.datetime(2015,12,3),var,"pe")
#postprocess_output(dt.datetime(2015,12,4),var,"pe")
#postprocess_output(dt.datetime(2015,12,5),var,"pe")
