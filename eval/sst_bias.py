import iris
from panMC import panMC
from bias_plots import bias_plots
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import numpy as np
import iris.plot as iplt
import subprocess
from iris.coord_categorisation import add_hour,add_month,add_day_of_year
from lagged_regressions import adjust_doyr
import datetime as dt
path = "/gws/nopw/j04/terramaris/panMC_um/{0}_*/postprocessed_outputs/sst/"
ref_path = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/"
sst12 = iris.cube.CubeList()
sstref = iris.cube.CubeList()

cx = iris.Constraint(longitude=lambda lon: 90<=lon<=155)
cy = iris.Constraint(latitude=lambda lat: -15<=lat<=15)

def load(year,MC):
  outname = "%s_%04d%02d_diurnal_sea_temperature.nc"%(MC,year,(year+1)%100)
  dcmean = iris.load(path.format(MC)+outname)[0]
  return dcmean


def load_ref(year,path,scratchpath=None,regen=False):
  outname = "ref_%04d%02d_surf_temperature.nc"%(year,(year+1)%100)
  if os.path.exists(path+"/monthlymean/"+outname) and not regen:
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

def daily_mean(year,MC):
  outname = "%s_%04d%02d_daily_sea_temperature.nc"%(MC,year,(year+1)%100)
  print(year)
  dmean = iris.load(path.format(MC)+outname)[0]
  import pdb;pdb.set_trace()
  return dmean

def load_ref_surf(year,path,level):
    data = iris.load(path+"%04d/%04d_00???_temperature_relax.nc"%(year,year))
    tmp = []
    times = []
    for cube in data:
      try:
        cube.coord("t").rename("time")
      except:
        1
      cube.coord("time").attributes={}
      cube.coord("latitude").attributes={}
      cube.coord("longitude").attributes={}
      cube.coord("time").long_name=None
      cube.coord("time").var_name=None
      n = cube.shape[0]
      for i in range(n):
        if cube[i].coord("time").points not in times:
          tmp.append(cube[i])
          times.append(cube[i].coord("time").points)
    iris.util.equalise_attributes(tmp)
    data=iris.cube.CubeList(tmp).merge_cube()
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    return data[:,level]

def calc_mse(year,MC,level,scratchpath,regen=False): 
  if os.path.exists(scratchpath+"%04d_%s_sst_mse.nc"%(year,MC)) and not regen:
     print("exists")
     mse = iris.load_cube(scratchpath+"%04d_%s_sst_mse.nc"%(year,MC))
     if mse.ndim==2:
       mse=mse = mse[:,level]
     return mse
  else:
    ct= iris.Constraint(time = lambda t: t.point.month in [12,1,2] and not (t.point.month,t.point.day)==(2,29))
    model =  daily_mean(year,MC)
    ref = load_ref_surf(year,ref_path,level).extract(ct)

    ref.coord("longitude").var_name=None
    model.coord("longitude").long_name=None
    model.coord("longitude").attributes={}
    ref.coord("longitude").attributes={}
    model.coord("longitude").guess_bounds()

    ref.coord("latitude").var_name=None
    model.coord("latitude").long_name=None
    model.coord("latitude").attributes={}
    ref.coord("latitude").attributes={}
    model.coord("latitude").guess_bounds()

    model.coord("time").convert_units(ref.coord("time").units)
    model.coord("time").bounds = None
    iris.util.equalise_attributes([model,ref])
    if MC=="MC2":
      ref = ref[:,54:-54,41:-41] 

    model = ((model - ref)**2).extract(cx&cy).collapsed(["longitude","latitude"],iris.analysis.MEAN)
    model.rename("mse_%s_sea_temperature"%MC)
    iris.save(model,scratchpath+"%04d_%s_sst_mse.nc"%(year,MC))
    return model

def combine_rmse(years,MC,level,scratchpath,regen=False):
 mse = iris.cube.CubeList()
 for year in years:
   yc = iris.coords.AuxCoord(year,long_name="year",units="years")
   mse.append(calc_mse(year,MC,level,scratchpath,regen))
   mse[-1].add_aux_coord(yc)
   adjust_doyr(mse[-1])
   mse[-1].remove_coord("time")
   iris.util.promote_aux_coord_to_dim_coord(mse[-1],"doyr")
 mse = mse.merge_cube()
 if mse.ndim == 2:
   rmse = mse.collapsed("year",iris.analysis.MEAN)**0.5
 else:
   rmse = mse**0.5
 return mse,rmse


def load_all(years12,years2,scratchpath=None,regen=False):
  sst2 = iris.cube.CubeList()
  sst12 = iris.cube.CubeList()
  sstref = iris.cube.CubeList()
  for year in years12:
    print(year)
    sstref.append(load_ref(year,ref_path,scratchpath,regen))
    sst12.append(load(year,"MC12"))
    sst12[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    sstref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  for year in years2:
    sst2.append(load(year,"MC2"))
    sst2[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  for cube in sst12+sst2:
    for x in cube.coords():
      x.var_name=None
  sst12=sst12.concatenate_cube()
  sst2=sst2.concatenate_cube()
  #sst12=sst12.concatenate_cube().aggregated_by("month",iris.analysis.MIN).collapsed("time",iris.analysis.MEAN)
  sstref=sstref.concatenate_cube()
  sstref.data = np.ma.masked_values(sstref.data,0)
  sst12.data.mask += sst12.data<=2
  sst2.data.mask += sst2.data<=2
  return sst12,sst2,sstref
 
def main(years12,years2,scratchpath,figname=None,regen=False):
  mse12, rmse12 = combine_rmse(years12,"MC12",5,scratchpath,regen=regen)
  mse2, rmse2 = combine_rmse(years2,"MC2",5,scratchpath,regen=regen)
  sst12,sst2,sstref = load_all(years12,years2,scratchpath,regen)
  if not figname is None:
    ct=iris.Constraint(time=lambda t:t.point.month in [12,1,2])
    sst12=sst12.collapsed("time",iris.analysis.MEAN)
    sst2=sst2.collapsed("time",iris.analysis.MEAN)
    #sst12=sst12.concatenate_cube().aggregated_by("month",iris.analysis.MIN).collapsed("time",iris.analysis.MEAN)
    sstref=sstref.extract(ct).collapsed("time",iris.analysis.MEAN)

    sstref.coord("longitude").var_name=None
    sstref.coord("longitude").attributes={}
    sstref.coord("latitude").var_name=None
    sstref.coord("latitude").attributes={}

    sst12.coord("longitude").long_name=None
    sst12.coord("longitude").attributes={}
    sst12.coord("longitude").guess_bounds()
    sst12.coord("latitude").long_name=None
    sst12.coord("latitude").attributes={}
    sst12.coord("latitude").guess_bounds()

    sst2.coord("longitude").long_name=None
    sst2.coord("longitude").attributes={}
    sst2.coord("longitude").guess_bounds()
    sst2.coord("latitude").long_name=None
    sst2.coord("latitude").attributes={}
    sst2.coord("latitude").guess_bounds()

    sst2 = sst2.regrid(sst12,iris.analysis.Linear(extrapolation_mode="mask"))
    #sst2_ = sst12.copy(data=np.ma.masked_array(np.ones(sst12.data.shape),mask=True))
    #sst2_.data[54:-54,41:-41]  = sst2.data
    #sst2 = sst2_

    fig=bias_plots([sst2[3],sst12[3],sstref[3]], ["MC2","MC12","Reference"] ,cmap1="magma",cmap2="bwr",projection=ccrs.PlateCarree(),minmaxN1rangeN2=(20,32,13,1.1,12),nmax=12,above=True,mark_inner=True)
    fig.set_figwidth(12)
    fig.set_figheight(7)
    fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
    fig.suptitle("Foundational Sea Surface Temperature")
    ax = fig.add_axes([0.07,0.07,0.5,0.2])
    ax.set_xticks([365-31,365-16,365+1,365+15,365+31,365+31+15,365+31+28])
    ax.set_xticklabels(["1-Dec","15-Dec","1-Jan","15-Jan","1-Feb","15-Feb","28-Feb"])
    iplt.plot(rmse12,label="MC12",axes=ax,c="k")
    iplt.plot(rmse2,label="MC2",axes=ax,c="b")
    [iplt.plot(mse12[i]**0.5,axes=ax,label=year,c="0.5",lw=1) for i,year in enumerate(years12)]
    ax.set_title("RMSE Error")
    #ax.set_xlim(365-31,365+31+28)
    fig.savefig(figname)
    import pdb;pdb.set_trace()

if __name__=="__main__":
  main(MC12_years,"/work/scratch-pw2/emmah/eval/","/home/users/emmah/eval_figs/sst_bias.png",regen=False)
