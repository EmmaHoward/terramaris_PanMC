#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python

import os
import iris
import datetime as dt
import numpy as np
from subprocess import check_call
import sys

t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1

t0 = dt.datetime(year,1,1)



print(7+(t1-t0).days)

path = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/u-cf309/%04d/%04d%02d%02dT0000Z/"%(year,t1.year,t1.month,t1.day)

outpath = "/work/scratch-nopw/emmah/postproc_N1280_fc/"
finalpath = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/%04d%02d_u-cf309/"%(year,(year+1)%100)


names = {"taux_in":"surface_downward_eastward_stress",
         "tauy_in":"surface_downward_northward_stress",
         "solar_in":"surface_net_downward_shortwave_flux",
         "nsolar_in":"surface_downward_heat_flux_in_sea_water_excludes_shortwave",
         "PminusE_in":"precipitation minus evaporation",
         "hmix":"ocean_mixed_layer_thickness",
         "S":"sea_water_salinity",
         "T":"sea_water_temperature",
         "fcorr_z":"heat_flux_correction_per_metre",
         "scorr":"virtual_salt_flux_correction_per_metre",
         "wT":"Turbulent_vertical_sea_water_temperature_flux",
         "wTnt":"Radiative_vertical_sea_water_temperature_flux",
         "wS":"Turbulent_vertical_salinity_flux",
         "rho":"sea_water_density",
         "cp":"specific_heat_capacity_of_sea_water"}

units = {"T":"degC",
         "S":"g/kg",
         "wT":"K.m/s",
         "wTnt":"K.m/s",
         "wS":"g/kg/s",
         "rho":"kg/m3",
         "cp":"J/kg/K",
         "scorr":"g/kg/s",
         "fcorr_z":"W/m3",
         "hmix":"m",
         "taux_in":"N/m2",
         "tauy_in":"N/m2",
         "solar_in":"W/m2",
         "nsolar_in":"W/m2",
         "PminusE_in":"mm/s"}


file_variables = {"sea_water_temperature":["T"],  # 1
         "sea_water_salinity":["S"],              # 2
         "heat_flux_correction":["fcorr_z"],      # 3
         "salinity_correction":["scorr"],         # 4
         "wT_turb":["wT"],                        # 5
         "wT_solar":["wTnt"],                     # 6
         "wS_turb":["wS"],                        # 7
         "sea_water_density":["rho"],             # 8
         "cp":["cp"],                             # 9
         "mixedlayerdepth":["hmix"],              #10
         "sea_surface":["taux_in","tauy_in","solar_in","nsolar_in","PminusE_in"]}


def postprocess_output(timestep,outname):
  variables = file_variables[outname] 
  data=iris.load(path+"KPPocean_%05d.nc"%timestep,variables)
  new = iris.cube.CubeList()
  for key in variables:
    cube=data.extract(key)[0]
    cube.rename(names[key])
    if cube.units=="unknown":
      cube.units = units[key] 
    else:
      cube.convert_units(units[key])
    cube.coord("longitude").units = "degrees"
    cube.coord("latitude").units = "degrees"
    [c for (c,i) in cube._dim_coords_and_dims if i==0][0].rename("time")
    if cube.ndim==4:
      [c for (c,i) in cube._dim_coords_and_dims if i==1][0].rename("depth")
    cube.coord("time").units="days since %04d-01-01T0000Z"%year
    cube.coord("time").convert_units("hours since %04d-11-01T0000Z"%year)
  t = t0+dt.timedelta(timestep-1)
  iris.save(data,outpath+"KPP/tma_N1280_KPPcoupled_kpp_%s_%04d%02d%02d.nc"%(outname,t.year,t.month,t.day),zlib=True)

def copy_remove(t,outname):
    # postprocess data from crun starting on date1 into file outname.nc
    # loop through times for file output
      outfile = "tma_N1280_KPPcoupled_kpp_%s_%04d%02d%02d.nc"%(outname,t.year,t.month,t.day)
      cmd = "mkdir -p %s/%s"%(finalpath,"kpp")
      check_call(cmd,shell=True)
      cmd = "cp %s/KPP/%s %s/kpp/%s"%(outpath,outfile,finalpath,outfile)
      check_call(cmd,shell=True)
      cmd = "rm %s/KPP/%s"%(outpath,outfile)
      print(cmd)
      check_call(cmd,shell=True)


job = int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
for timestep in range(1+(t1-t0).days,7+(t1-t0).days):
#  for job in range(10):
    var = list(file_variables.keys())[job]
    print(var,timestep)
    postprocess_output(timestep,var)
    copy_remove(t0+dt.timedelta(timestep-1),var)
