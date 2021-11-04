#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#  SBATCH -p test
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH --array=[1-48]
#SBATCH -o /home/users/emmah/log/postproc_2km/pp2_%a.o
#SBATCH -e /home/users/emmah/log/postproc_2km/pp2_%a.e 
#SBATCH -t 04:00:00
from subprocess import check_call
import sys
import iris
import datetime as dt
import os

t0 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))
if t0.month >6:
  year = t0.year
else:
  year = t0.year - 1

path = "/work/scratch-nopw/emmah/pepf/%04d%02d%02dT0000Z/"%(t0.year,t0.month,t0.day)
check_call("mkdir -p %s/sep1"%path,shell=True)
os.chdir(path)

c_mean = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method  in ["mean","sum"]))
c_max = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="maximum"))
c_min = iris.Constraint(cube_func=lambda cube: len(cube.cell_methods)>0 and (cube.cell_methods[0].method=="minimum"))
c_inst = iris.Constraint(cube_func=lambda cube: not (cube.cell_methods) or (cube.cell_methods[0].method=="point"))


c3 = iris.Constraint(cube_func=lambda cube: cube.ndim==3)
c4 = iris.Constraint(cube_func=lambda cube: cube.ndim==4)


file_variables_pepf = {
         "cloud_fraction"       :([c_inst],["cloud_area_fraction_in_atmosphere_layer",                               # 1
                                         "cloud_volume_fraction_in_atmosphere_layer",
                                         "liquid_cloud_volume_fraction_in_atmosphere_layer",
                                         "ice_cloud_volume_fraction_in_atmosphere_layer"]),
#                                         "ice_aggregate_fraction"]),

         "radar_reflectivity"   :([c_inst],["radar_reflectivity_due_to_graupel_alone",                               # 2
                                          "radar_reflectivity_due_to_ice_aggregates_alone",
#                                          "radar_reflectivity_due_to_ice_crystals_alone",
                                          "radar_reflectivity_due_to_rain_alone",
                                          "radar_reflectivity_due_to_cloud_alone"]),
         "density"              :([c_inst],["air_density"]),                                                         # 3
         "potential_temperature":([c_inst],["air_potential_temperature"]),                                           # 4
         "pressure_rho_grid"    :([c_inst],["m01s00i407"]),                                                          # 5
         "pressure_theta_grid"  :([c_inst],["m01s00i408"]),                                                          # 6
         "mixing_ratios"        :([c_inst],["mass_fraction_of_cloud_ice_in_air",                                     # 7
                                                 "mass_fraction_of_cloud_liquid_water_in_air",
                                                 "mass_fraction_of_graupel_in_air",
                                                 "mass_fraction_of_rain_in_air"]),
         "specific_humidity"    :([c_inst],["specific_humidity"]),                                                   # 8
         "eastward_wind"        :([c_inst],["x_wind"]),                                                         # 9
         "northward_wind"       :([c_inst],["y_wind"]),
         "upward_air_velocity"  :([c_inst],["upward_air_velocity"])}


def separate(t,stream):
  data=iris.load(path+"tm*a_%s%04d%02d%02d_%02d??"%(stream,t.year,t.month,t.day,t.hour))
  for key in ["cloud_fraction","radar_reflectivity","density","potential_temperature","pressure_rho_grid","pressure_theta_grid","mixing_ratios","specific_humidity","eastward_wind","northward_wind","upward_air_velocity"]:
    c_agg, vars = file_variables_pepf[key]
    subset = data.extract(c4).extract(vars).extract(c_agg)
    print(subset)
    assert len(vars)==len(subset)
    iris.save(subset,"sep1/%04d%02d%02d%02d_%s_%s.nc"%(t.year,t.month,t.day,t.hour,key,stream))  
  subset = data.extract(c3)
  iris.save(subset,"sep1/%04d%02d%02d%02d_%s_%s.nc"%(t.year,t.month,t.day,t.hour,"single_level",stream))
  check_call("rm %s/tm*a_%s%04d%02d%02d_%02d??"%(path,stream,t.year,t.month,t.day,t.hour),shell=True)

def separate_fill(t,stream):
  for key in ["cloud_fraction","radar_reflectivity","density","potential_temperature","pressure_rho_grid","pressure_theta_grid","mixing_ratios","specific_humidity","eastward_wind","northward_wind","upward_air_velocity"]:
    if not os.path.exists("sep1/%04d%02d%02d%02d_%s_%s.nc"%(t.year,t.month,t.day,t.hour,key,stream)):
      data=iris.load(path+"tm*a_%s%04d%02d%02d_%02d??"%(stream,t.year,t.month,t.day,t.hour))
      c_agg, vars = file_variables_pepf[key]
      subset = data.extract(c4).extract(vars).extract(c_agg)
      print(subset)
      iris.save(subset,"sep1/%04d%02d%02d%02d_%s_%s.nc"%(t.year,t.month,t.day,t.hour,key,stream))  
  if not os.path.exists("sep1/%04d%02d%02d%02d_%s_%s.nc"%(t.year,t.month,t.day,t.hour,"single_level",stream)):
    data=iris.load(path+"tm*a_%s%04d%02d%02d_%02d??"%(stream,t.year,t.month,t.day,t.hour))
    subset = data.extract(c3)
    iris.save(subset,"sep1/%04d%02d%02d%02d_%s_%s.nc"%(t.year,t.month,t.day,t.hour,"single_level",stream))


  
def untar(i):
  cycle="%04d%02d%02dT0000Z"%(t0.year,t0.month,t0.day)
  gwspath="/gws/nopw/j04/terramaris/emmah/coupled_2km/u-cc339/2015/"
  scratchpath="/work/scratch-nopw/emmah/pepf"
  check_call("mkdir -p %s/%s/sep1"%(scratchpath,cycle),shell=True)
  check_call("cp %s/%s/tm2a_pe_%03d.tar %s/%s"%(gwspath,cycle,i,scratchpath,cycle),shell=True)
  check_call("cp %s/%s/tm2a_pf_%03d.tar %s/%s"%(gwspath,cycle,i,scratchpath,cycle),shell=True)
  check_call("tar -xvf tm2a_pe_%03d.tar"%i,shell=True)
  check_call("tar -xvf tm2a_pf_%03d.tar"%i,shell=True)
  check_call("rm %s/%s/tm2a_pe_%03d.tar"%(scratchpath,cycle,i),shell=True)
  check_call("rm %s/%s/tm2a_pf_%03d.tar"%(scratchpath,cycle,i),shell=True)

job =  int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
untar(job+1)
for i in range(3):
  t = t0+dt.timedelta((job*3+i)/24)
  separate(t,"pe")
  separate(t,"pf")
