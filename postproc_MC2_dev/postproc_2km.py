#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python

import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call
from iris.aux_factory import HybridHeightFactory
import sys

print(sys.argv[1])
i_start = int(sys.argv[2])
i_end = int(sys.argv[3])

t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1


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
if t1.day ==1 and i_start in [0,1]:
  reinit_step_pp = {"pa":4,"pb":4,"pc":4,"pd":6} # input, in hours
else:
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
        return fac
 

def postprocess_output(date1,outname):
    stream,cell_methods,variables = file_variables[outname] 
    # postprocess data from crun starting on date1 into file outname.nc
    # loop through times for file output
    if stream == "pa":
       reinit = reinit_step_pa[outname]
    else:
       reinit = reinit_step_nc[stream]
    for i in range(i_start*24,i_end*24,reinit):
      date2 = date1 + dt.timedelta(i/24)
      if date1.day==1:
        path = "/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/%04d/%04d%02d%02dT0000Z/"%(year,date2.year,date2.month,date2.day)
      else:
        path = "/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/%04d/%04d%02d%02dT0000Z/"%(year,date1.year,date1.month,date1.day)
      print(date2)
      # times for file input
      times = [date2+ dt.timedelta(j/24) for j in np.arange(0,reinit,reinit_step_pp[stream])]
      files = [path+"/tm2a_%s%04d%02d%02d_%02d*"%(stream,t.year,t.month,t.day,t.hour)\
             for t in times]
      data = iris.load(files,callback=callback)#.extract(cx&cy&cz)
      for cube in data:
        if cube.name() in units.keys():
            # correct units
            cube.units = units[cube.name()]
        if cube.name() in names.keys():
            # correct variable names
            cube.rename(names[cube.name()])
      if stream=="pd" and cube.name()!="Heaviside function on pressure levels":
         # extract heaviside function for masking data below orography
         H = data.extract("Heaviside function on pressure levels")[0].copy()
      # extract variables for saving 
      if outname in ["w_T_rho","w_q_rho"]:
        fac = get_wxrho_factor(data)
      data = data.extract(variables).extract(cell_methods)
      print(data)
      for cube in data:
        # convert units to hour since simulation started
        cube.coord("time").convert_units("hours since 2015-11-01 00:00:00")
        # convert to float
        cube.data = cube.data.astype("float32")
        for c in cube.coords():
          if c.points.dtype in [int,"int32"]:
            c.points = c.points.astype("int32")
            if c.has_bounds():
              c.bounds = c.bounds.astype("int32")
          else:
            c.points = c.points.astype("float32")
            if c.has_bounds():
              c.bounds = c.bounds.astype("float32")         
          if c.name() in ["latitude","longitude"]:# and stream != "pc":
            c.points = np.round(c.points * 100)/100
      if stream in ["pc"]:
        # for pc stream: correct wxrho terms and coarsen to N1280 grid
        new = iris.cube.CubeList()
        for cube in data:
          add_altitude_coord(cube)
          if outname in ["w_T_rho","w_q_rho"]:
            z_r = cube.coord("altitude").points
            earth_radius = 6371229.0
            cube.data = cube.data*fac*(earth_radius**2)/((z_r + earth_radius)**2)
          new.append(interp(cube))
        data=new
      data=data.concatenate()
      data=data.merge()
      # check for the right number of variables. "range" has 2x "variables" due to max and mins, except wind gust
      if outname == "range": 
        assert len(data) == len(variables)*2-1
      else:
        assert len(data) == len(variables)
      # save to netcdf
      if reinit==24:
        outfile = "tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"%(stream,outname,date2.year,date2.month,date2.day)
      else:
        outfile = "tma_2km_KPPcoupled_%s_%s_%04d%02d%02d_%02d.nc"%(stream,outname,date2.year,date2.month,date2.day,date2.hour)
      iris.save(data,outpath+"%s/%s"%(stream,outfile),zlib=True)

def collate_pc(date1,outname):
   stream,cell_methods,variables = file_variables[outname]
   reinit = reinit_step_nc[stream]
   for j in range(i_start,i_end):
     outfiles = []
     date3 = date1+dt.timedelta(j)
     for i in range(j*24,(j+1)*24,reinit):
       date2 = date1 + dt.timedelta(i/24)
       outfiles.append(outpath+"/%s/tma_2km_KPPcoupled_%s_%s_%04d%02d%02d_%02d.nc"%(stream,stream,outname,date2.year,date2.month,date2.day,date2.hour))
     data=iris.load(outfiles,variables)
     data=data.concatenate_cube()
     iris.save(data,outpath+"/%s/tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"%(stream,stream,outname,date3.year,date3.month,date3.day),zlib=True)
     #for outfile in outfiles:
     #  check_call("rm %s"%outfile,shell=True)

def callback(cube, field, filename):
  cube.remove_coord("forecast_period")
  cube.remove_coord("forecast_reference_time")

def copy_ncks_remove(date1,outname):
    # postprocess data from crun starting on date1 into file outname.nc
    stream,cell_methods,variables = file_variables[outname] 
    # loop through times for file output
    if stream == "pa":
       reinit = reinit_step_pa[outname]
    else:
       reinit = reinit_step_nc[stream]
    for i in range(i_start*24,i_end*24,reinit):
      date2 = date1 + dt.timedelta(i/24)
      if reinit == 24:
        outfile = "tma_2km_KPPcoupled_%s_%s_%04d%02d%02d.nc"%(stream,outname,date2.year,date2.month,date2.day)
      else:
        outfile = "tma_2km_KPPcoupled_%s_%s_%04d%02d%02d_%02d.nc"%(stream,outname,date2.year,date2.month,date2.day,date2.hour)
      print(outfile)
      os.chdir(outpath+stream)
      cmd = "ncks -L 1 %s ncks/%s"%(outfile,outfile)
      if 1:#not os.path.exists("%s/%s/%s"%(finalpath,stream,outfile)):
        if 1:# not os.path.exists("ncks/%s"%(outfile)):
          check_call(cmd,shell=True)
        cmd = "cp ncks/%s %s/%s/%s"%(outfile,finalpath,stream,outfile)
        check_call(cmd,shell=True)
      cmd = "rm %s ncks/%s"%(outfile,outfile)
      check_call(cmd,shell=True)


job =  int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
var = list(file_variables.keys())[job]
print(var)
postprocess_output(t1,var)
if file_variables[var][0]=="pc":
    collate_pc(t1,var)
    reinit_step_nc["pc"]=24
copy_ncks_remove(t1,var)











