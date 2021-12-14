import iris.plot as iplt
import iris
from panMC import panMC
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import subprocess
from iris.coord_categorisation import add_hour,add_month
MC12_years = [2003,2014,2015,2016,2017]
MC12_path = "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/monthly_diurnal_sst/"
ref_path = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/"
scratchpath = "/work/scratch-nopw/emmah/"
sst12 = iris.cube.CubeList()
sstref = iris.cube.CubeList()
def load(year,MC,path):
  outname = "%s_%04d%02d_sea_temperature.nc"%(MC,year,(year+1)%100)
  print(year)
  if os.path.exists(path+outname):
    print("exists")
    dcmean = iris.load(path+outname)[0]
  else:
    data=panMC(year,MC,"sea_water_temperature").load_iris()[0]
    import pdb;pdb.set_trace()
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    data.coord("time").bounds = np.array([np.round(data.coord("time").points)-1,np.round(data.coord("time").points)]).T
    data.coord("time").points = np.round(data.coord("time").points)-0.5
    add_hour(data,"time","hour")
    add_month(data,"time","month")
    dcmean = data[:,:6].aggregated_by(["hour","month"],iris.analysis.MEAN)
    iris.save(dcmean,scratchpath+outname,zlib=True)
    subprocess.check_call("cp %s%s %s%s"%(scratchpath,outname,path,outname),shell=True)
    subprocess.check_call("rm %s%s"%(scratchpath,outname),shell=True)
  return dcmean


def load_ref(year,path):
  outname = "ref_%04d%02d_surf_temperature.nc"%(year,(year+1)%100)
  if os.path.exists(path+"/monthlymean/"+outname):
    print("exists")
    mean = iris.load(path+"/monthlymean/"+outname)[0]
  else:
    data = iris.load(path+"%04d/%04d_00???_temperature_relax.nc"%(year,year))
    data = iris.cube.CubeList([cube[:-1] for cube in data if cube.shape[0]==7]).concatenate_cube()
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    data.coord("t").rename("time")
    add_month(data,"time","month")  
    mean = data[:,:6].aggregated_by("month",iris.analysis.MEAN)
    iris.save(mean,scratchpath+outname,zlib=True)
    subprocess.check_call("cp %s%s %s/monthlymean/%s"%(scratchpath,outname,path,outname),shell=True)
    subprocess.check_call("rm %s%s"%(scratchpath,outname),shell=True)
  return mean 

ct=iris.Constraint(time=lambda t:t.point.month in [12,1,2])
for year in MC12_years:
  print(year)
  sstref.append(load_ref(year,ref_path).extract(ct).collapsed("time",iris.analysis.MEAN))
  sst12.append(load(year,"MC12",MC12_path).extract(ct).collapsed("time",iris.analysis.MEAN))
  sst12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  sstref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  sstref[-1].data = np.ma.masked_array(sstref[-1].data,mask = sstref[-1].data < 2)



for cube in sst12:
   for x in cube.coords():
      x.var_name=None
  
sst12=sst12.merge_cube().collapsed("longitude",iris.analysis.MEAN)
#sst12=sst12.concatenate_cube().aggregated_by("month",iris.analysis.MIN).collapsed("time",iris.analysis.MEAN)
sstref=sstref.merge_cube().collapsed("longitude",iris.analysis.MEAN)


fig=plt.figure(figsize=(5,7))
plt.subplot(311)
plt.title("MC2")
plt.xlim(-20,20)
plt.plot([-15,-15,15,15,-15],[24,31,31,24,24],"k",lw=1)
plt.ylabel("SST (C)")
plt.ylim(25,30)
plt.subplot(312)
plt.plot([-15,-15,15,15,-15],[24,31,31,24,24],"k",lw=1)
plt.title("MC12")
plt.ylabel("SST (C)")
plt.ylim(25,30)
plt.xlim(-20,20)
for i,year in enumerate(MC12_years):
  iplt.plot(sst12[i,3])#,label="%04d-%02d"%(year,(year+1)%100))

plt.subplot(313)
plt.plot([-15,-15,15,15,-15],[24,31,31,24,24],"k",lw=1)
plt.xlim(-20,20)
plt.ylim(25,30)
plt.title("Ostia SST")
plt.ylabel("SST (C)")
for i,year in enumerate(MC12_years):
  iplt.plot(sstref[i,3],label="%04d-%02d"%(year,(year+1)%100))

fig.subplots_adjust(right=0.8,hspace=0.3)
fig.legend(loc="right")

plt.show()
