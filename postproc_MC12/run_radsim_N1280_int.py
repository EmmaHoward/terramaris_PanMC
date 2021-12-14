#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python

import iris
import datetime as dt
import numpy as np
from subprocess import check_call
import os
import sys
t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))
if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1


c_inst = iris.Constraint(cube_func=lambda cube: not cube.cell_methods)

# simulation start
t0 = dt.datetime(year,11, 1)
# crun start
path = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/u-cf309/%04d/%04d%02d%02dT0000Z/"%(t0.year,t1.year,t1.month,t1.day)

#keys for multilevel fields needed for radsim
keys_a = ["air_potential_temperature","specific_humidity","m01s00i408","cloud_area_fraction_in_atmosphere_layer","mass_fraction_of_cloud_ice_in_air","mass_fraction_of_cloud_liquid_water_in_air"]

#keys for single level fields needed for radsim
keys_b = ["surface_temperature","m01s03i225","m01s03i226","air_temperature","relative_humidity","surface_air_pressure"]

orogpath = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"

# input re-initialisation timestep in hours
reinit_step_pp = {"pa":4,"pb":24,"pc":6,"pd":12,"pg":4}

def prep_data_pg(date):
  # builds an input pp file for one timestep of radsim
  ct = iris.Constraint(time = lambda t: (t.point.year==date.year and t.point.month==date.month and t.point.day==date.day and t.point.hour == date.hour and t.point.minute==date.minute))
  #orog = iris.load(path+"tma_pb%04d%02d%02d_00.pp"%(t1.year,t1.month,t1.day),["surface_altitude"])
  orog = iris.load(orogpath)
  pg_hr = ((date.hour-(date.minute==0)) // reinit_step_pp["pg"])*reinit_step_pp["pg"]
  # read data
  if date.hour == 0 and date.minute == 0:
    # hour 0 is in file labelled with previous day
    date1=date+dt.timedelta(-1)
    pg = iris.load(path+"tma_pg%04d%02d%02d_%02d.pp"%(date1.year,date1.month,date1.day,20),ct).extract(keys_a+keys_b)
  else:
    pg = iris.load(path+"tma_pg%04d%02d%02d_%02d.pp"%(date.year,date.month,date.day,pg_hr),ct).extract(keys_a+keys_b)
  assert(len(pg)==len(keys_a+keys_b))
  # seaice is zero everywhere and not outputted - make a dummy variable for it
  seaice = pg.extract("surface_temperature")[0].copy(data = np.zeros(pg.extract("surface_temperature")[0].shape))
  seaice.rename("fractional sea ice")
  seaice.units=(1)
  clf = pg.extract("cloud_area_fraction_in_atmosphere_layer")[0]
  clf.data = np.maximum(clf.data,0)
  clf.data = np.minimum(clf.data,1)
  seaice.attributes["STASH"] = iris.fileformats.pp.STASH(1,0,31)
  for cube in pg+orog:
    cube.remove_coord("forecast_reference_time")
    cube.remove_coord("forecast_period")
  # orography is time invariant
  #orog[0].remove_coord("time")
  #orog[0].add_aux_coord(pb[0].coord("time"))
  # save to pp file
  print("/work/scratch-nopw/emmah/radsim/radsim_%04d%02d%02d_%02d%02d.pp"%(date.year,date.month,date.day,date.hour,date.minute))
  iris.save(pg+orog+[seaice],"/work/scratch-nopw/emmah/radsim/radsim_%04d%02d%02d_%02d%02d.pp"%(date.year,date.month,date.day,date.hour,date.minute))


def prep_data(date):
  # builds an input pp file for one timestep of radsim
  ct = iris.Constraint(time = lambda t: (t.point.year==date.year and t.point.month==date.month and t.point.day==date.day and t.point.hour == date.hour))
  #orog = iris.load(path+"tma_pb%04d%02d%02d_00.pp"%(t1.year,t1.month,t1.day),["surface_altitude"])
  orog = iris.load(orogpath)
  pa_hr = ((date.hour-1) // reinit_step_pp["pa"])*reinit_step_pp["pa"]
  # read data
  if date.hour == 0:
    # hour 0 is in file labelled with previous day
    date1=date+dt.timedelta(-1)
    pa = iris.load(path+"tma_pa%04d%02d%02d_%02d.pp"%(date1.year,date1.month,date1.day,20),ct).extract(keys_a)
    pb = iris.load(path+"tma_pb%04d%02d%02d_%02d.pp"%(date1.year,date1.month,date1.day,0),ct&c_inst).extract(keys_b) 
  else:
    pa = iris.load(path+"tma_pa%04d%02d%02d_%02d.pp"%(date.year,date.month,date.day,pa_hr),ct).extract(keys_a)
    pb = iris.load(path+"tma_pb%04d%02d%02d_%02d.pp"%(date.year,date.month,date.day,0),ct&c_inst).extract(keys_b) 
  assert(len(pa)==len(keys_a))
  assert(len(pb)==len(keys_b))
  # seaice is zero everywhere and not outputted - make a dummy variable for it
  seaice = pb.extract("surface_temperature")[0].copy(data = np.zeros(pb.extract("surface_temperature")[0].shape))
  seaice.rename("fractional sea ice")
  seaice.units=(1)
  seaice.attributes["STASH"] = iris.fileformats.pp.STASH(1,0,31)
  for cube in pa+pb+orog:
    cube.remove_coord("forecast_reference_time")
    cube.remove_coord("forecast_period")
  # orography is time invariant
  #orog[0].remove_coord("time")
  #orog[0].add_aux_coord(pb[0].coord("time"))
  # save to pp file
  print("/work/scratch-nopw/emmah/radsim/radsim_%04d%02d%02d_%02d%02d.pp"%(date.year,date.month,date.day,date.hour,date.minute))
  iris.save(pa+pb+orog+[seaice],"/work/scratch-nopw/emmah/radsim/radsim_%04d%02d%02d_%02d%02d.pp"%(date.year,date.month,date.day,date.hour,date.minute))


def build_namelist(date):
  # copies RadSim namelist with appropriate date
  # Read in the file
  with open('/home/users/emmah/radsim/radsim_cfg_base.nl', 'r') as file :
    filedata = file.read()
  # Replace the target string
  filedata = filedata.replace('YYYYMMDD_HHmm', '%04d%02d%02d_%02d%02d'%(date.year,date.month,date.day,date.hour,date.minute))
  # Write the file out again
  with open("/work/scratch-nopw/emmah/radsim/radsim_cfg_%04d%02d%02d_%02d%02d.nl"%(date.year,date.month,date.day,date.hour,date.minute), 'w') as file:
    file.write(filedata)
  file.close()
  print("written to /work/scratch-nopw/emmah/radsim/radsim_cfg_%04d%02d%02d_%02d%02d.nl"%(date.year,date.month,date.day,date.hour,date.minute))

def run_radsim(date):
  # run radiance simulator
  namelist = "/work/scratch-nopw/emmah/radsim/radsim_cfg_%04d%02d%02d_%02d%02d.nl"%(date.year,date.month,date.day,date.hour,date.minute)
  cmd = "/home/users/emmah/radsim/radsim-2.2/bin/radsim.exe %s"%namelist
  check_call(cmd, shell=True)

def pp_radsim(date):
  # postprocess output from radiance simulator
  # load data
  data=iris.load("/work/scratch-nopw/emmah/radsim/radsim-himawari_8_ahi-%04d%02d%02d%02d%02d.nc"%(date.year,date.month,date.day,date.hour,date.minute))
  for key in ["qcrttov","qcinfo","qcflags"]:
   # check quality control flags, write warning file if any flags were raised
   if (data.extract("qcrttov")[0].data.data!=0).any():
     print("qc flag raised")
     f = open("/work/scratch-nopw/emmah/radsim/%04d%02d%02d%02d%02d_qcflag.txt"%(date.year,date.month,date.day,date.hour,date.minute),"w")
     f.write(key,np.where(data.extract("qcrttov")[0].data.data!=0))
     f.close()
  # extract brightness temperatures
  bt = data.extract("bt")[0].data
  # load orography to obtain terramaris grid
  orog = iris.load(orogpath)[0]
  # determine number of channels
  nc = bt.shape[0]
  bt = bt.reshape(nc,orog.shape[0],orog.shape[1])
  # create iris cube
  bt = iris.cube.Cube(bt,standard_name="brightness_temperature",units="K")
  bt.add_dim_coord(orog.coord("latitude"),1)
  bt.add_dim_coord(orog.coord("longitude"),2)
  # create sensor band, wavenumber and frequency coordinates
  bands = iris.coords.DimCoord(np.arange(7,17).astype(float),standard_name="sensor_band_identifier",units=1)
  wn = iris.coords.AuxCoord(data.extract("bt")[0].attributes["wavenumbers"],standard_name="sensor_band_central_radiation_wavenumber",units="cm-1")
  freq = iris.coords.AuxCoord(1/data.extract("bt")[0].attributes["wavenumbers"],standard_name="sensor_band_central_radiation_wavelength",units="cm")
  freq.convert_units("micron")
  bt.add_dim_coord(bands,0)
  bt.add_aux_coord(wn,0)
  bt.add_aux_coord(freq,0)
  # check time
  assert (np.array([date.year,date.month,date.day,date.hour,date.minute])  == data.extract("bt")[0].attributes["validity_time"]).all()
  dateC = iris.coords.DimCoord((date - t0).total_seconds()/3600,standard_name="time",units="hours since %04d-%02d-%02d %02d:00:00"%(t0.year,t0.month,t0.day,t0.hour))
  bt.add_aux_coord(dateC)
  # copy metadata
  for key in ["instrument","platform","satid"]:
    bt.attributes[key] = data.extract("bt")[0].attributes[key]
  # save to netcdf file
  iris.save(bt,"/work/scratch-nopw/emmah/radsim/tma_radsim-himawari_8_ahi-%04d%02d%02d_%02d%02d.nc"%(date.year,date.month,date.day,date.hour,date.minute),zlib=True)

def main(date):
  print(date)
  # build namelist for radsim with correct date
  build_namelist(date)
  # extract relevant model output
  if date.minute==0:
    prep_data(date)
  else:
    prep_data_pg(date)
  # run radiance simulator
  run_radsim(date)
  # postprocess output
  pp_radsim(date)
  # remove intermediate files
  if not os.path.exists("/work/scratch-nopw/emmah/radsim/%04d%02d%02d%02d%02d_qcflag.txt"%(date.year,date.month,date.day,date.hour,date.minute)):
    check_call("rm /work/scratch-nopw/emmah/radsim/radsim-himawari_8_ahi-{year:04d}{month:02d}{day:02d}{hour:02d}{mins:02d}.nc \
                   /work/scratch-nopw/emmah/radsim/radsim_cfg_{year:04d}{month:02d}{day:02d}_{hour:02d}{mins:02d}.nl \
                   /work/scratch-nopw/emmah/radsim/radsim_{year:04d}{month:02d}{day:02d}_{hour:02d}{mins:02d}.pp".format(year=date.year,month=date.month,day=date.day,hour=date.hour,mins=date.minute),shell=True)

"""
job = int(os.environ["SLURM_ARRAY_TASK_ID"])-1
#for i in range(4,145):
for i in range(job*3,job*3+3,1):
  date = (t1+dt.timedelta(i/24))
  for minute in [0,21,39]:
    if minute==0 and i==0:
      date1 = date+dt.timedelta(6)
    else:
      date1 =date.replace(minute=minute)
    if not (os.path.exists("/work/scratch-nopw/emmah/radsim/tma_radsim-himawari_8_ahi-%04d%02d%02d_%02d%02d.nc"%(date1.year,date1.month,date1.day,date1.hour,date1.minute))):
      print(date1)
      main(date1)
  
"""
# uncomment this and comment out above for loop to identify and interactively run timesteps that were missed by lotus 
for i in range(0,24*6):
  date = (t1+dt.timedelta(i/24))
  for minute in [0,21,39]:
    if minute==0 and i==0:
      date1 = date+dt.timedelta(6)
    else:
      date1 =date.replace(minute=minute)
    if not (os.path.exists("/work/scratch-nopw/emmah/radsim/tma_radsim-himawari_8_ahi-%04d%02d%02d_%02d%02d.nc"%(date1.year,date1.month,date1.day,date1.hour,date1.minute))):
      print(date1)
      main(date1)
