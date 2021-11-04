#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#SBATCH -p short-serial-4hr
#  SBATCH -p high-mem
#SBATCH -A short4hr
#SBATCH --array=[1,6,7,8,9,10]
#SBATCH -o /home/users/emmah/log/postproc_2km/ncks1_%a.o
#SBATCH -e /home/users/emmah/log/postproc_2km/ncks1_%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=72000

import os
from subprocess import check_call

file_variables = {
         "radar_reflectivity"   :"pa",   # 1
         "density"              :"pa",   # 2
         "potential_temperature":"pa",   # 3
         "pressure_rho_grid"    :"pa",   # 4
         "pressure_theta_grid"  :"pa",   # 5
         "cloud_fraction"       :"pa",   # 6
         "qcf"                  :"pa",   # 7
         "qcl"                  :"pa",   # 8
         "rain"                 :"pa",   # 9
         "graupel"              :"pa",   #10
         "specific_humidity"    :"pa",   #11
         "upward_air_velocity"  :"pa",   #12
         "eastward_wind"        :"pa",   #13
         "northward_wind"       :"pa",   #14
#
         "atmos":                "pb",   #15
         "integral_inst":        "pb",   #16
         "integral_agg":         "pb",   #17
         "rainfall":             "pb",   #18
         "surf_inst":            "pb",   #19
         "surf_agg":             "pb",   #20
         "range":                "pb",   #21
#                                         
         "dT_shortwave"          :"pc",  #22
         "dT_longwave"           :"pc",  #23
         "dT_advection"          :"pc",  #24
         "dq_advection"          :"pc",  #25
         "w_T_rho"               :"pc",  #26
         "w_q_rho"               :"pc",  #27
         "dT_total"              :"pc",  #28
         "dq_total"              :"pc",  #29
         "temperature_pc"        :"pc",  #30
         "specific_humidity_pc"  :"pc",  #31
         "upward_air_velocity_pc":"pc",  #32 
#
         "Heavyside"             :"pd",  #33
         "temperature_pd"        :"pd",  #34
         "geopotential_height"   :"pd",  #35
         "omega"                 :"pd",  #36
         "relative_humidity"     :"pd",  #37
         "specific_humidity_pd"  :"pd",  #38
         "theta_w"               :"pd",  #39
         "eastward_wind_pd"      :"pd",  #40
         "div_vort"              :"pd",  #41
         "northward_wind_pd"     :"pd",  #42
#         "potential_vorticity"   :("pd",[c_inst],["potential_vorticity"]),                                                     #
         "daily_mean"            :"pd"   #43
         
         }     


job = int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
var = list(file_variables.keys())[job]
print(var)
stream = file_variables[var]
outpath = "/work/scratch-nopw/emmah/postproc_2km_fc/"
os.chdir(outpath+stream)
#files = [f for f in os.listdir(outpath+stream) if var in f]

f = "tma_2km_KPPcoupled_pa_%s_20151201.nc"%var
cmd = "ncks -L 1 %s ncks/%s"%(f,f)
print(cmd)
check_call(cmd,shell=True)


