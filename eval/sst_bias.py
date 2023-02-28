#
# plot sst biases (Figure 2)
#

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

# domain bounds
cx = iris.Constraint(longitude=lambda lon: 90.5<=lon<=154.5)
cy = iris.Constraint(latitude=lambda lat: -15<=lat<=15)

def load(year,MC):
  # load sst data (diurnal mean)
  outname = "%s_%04d%02d_diurnal_sea_temperature.nc"%(MC,year,(year+1)%100)
  dcmean = iris.load(path.format(MC)+outname)[0]
  return dcmean


def load_ref(year,path,scratchpath=None,regen=False):
  # load reference sst data (return monthly mean)
  outname = "ref_%04d%02d_surf_temperature.nc"%(year,(year+1)%100)
  if os.path.exists(path+"/monthlymean/"+outname) and not regen:
    print("exists")
    mean = iris.load(path+"/monthlymean/"+outname)[0]
  else:
    data = iris.load(path+"%04d/%04d_00???_temperature_relax.nc"%(year,year))
    times,tmp = [],iris.cube.CubeList()
    for cube in data:
      if "t" in [c.name() for c in cube.coords()]:
        time = cube.coord("t").units.num2date(cube.coord("t").points)
      else:
        time = cube.coord("time").units.num2date(cube.coord("time").points)
      for i,t in enumerate(time):
        if t not in times:
           tmp.append(cube[i])
           times.append(t)
    data = tmp.merge_cube()
    [c for (c,i) in data._dim_coords_and_dims if i==1][0].rename("depth")
    data.coord("t").rename("time")
    add_month(data,"time","month")  
    mean = data[:,:6].aggregated_by("month",iris.analysis.MEAN)
    iris.save(mean,scratchpath+outname,zlib=True)
    subprocess.check_call("cp %s%s %s/monthlymean/%s"%(scratchpath,outname,path,outname),shell=True)
    subprocess.check_call("rm %s%s"%(scratchpath,outname),shell=True)
  return mean 

def daily_mean(year,MC):
  # load sst data (daily mean)
  outname = "%s_%04d%02d_daily_sea_temperature.nc"%(MC,year,(year+1)%100)
  print(year)
  dmean = iris.load(path.format(MC)+outname)[0]
  return dmean

def load_ref_surf(year,path,level):
    # load reference sst data (return single model level, daily data)
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


"""
old: interannual rmse

def calc_mse(year,MC,level,scratchpath,regen=False): '
  # compute mean square error for a single year
  if os.path.exists(scratchpath+"%04d_%s_sst_mse.nc"%(year,MC)) and not regen:
     # load cached data
     print("exists")
     mse = iris.load_cube(scratchpath+"%04d_%s_sst_mse.nc"%(year,MC))
     if mse.ndim==2:
       mse=mse = mse[:,level]
     return mse
  else:
    # load daily data
    ct= iris.Constraint(time = lambda t: t.point.month in [12,1,2] and not (t.point.month,t.point.day)==(2,29))
    model =  daily_mean(year,MC)[:,level]
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
      # restrict to MC2 domain
      ref = ref[:,54:-54,41:-41] 
    # compute mse
    model = ((model - ref)**2).extract(cx&cy).collapsed(["longitude","latitude"],iris.analysis.MEAN)
    model.rename("mse_%s_sea_temperature"%MC)
    # cache data
    iris.save(model,scratchpath+"%04d_%s_sst_mse.nc"%(year,MC))
    return model


def combine_rmse(years,MC,level,scratchpath,regen=False):
 # compute overall rmse from annual mses
 mse = iris.cube.CubeList()
 for year in years:
   yc = iris.coords.AuxCoord(year,long_name="year",units="years")
   mse.append(calc_mse(year,MC,level,scratchpath,regen))
   mse[-1].add_aux_coord(yc)
   adjust_doyr(mse[-1])
   mse[-1].remove_coord("time")
   if "depth" in [c.name() for c in mse[-1].coords()]:
     mse[-1].remove_coord("depth")
   iris.util.promote_aux_coord_to_dim_coord(mse[-1],"doyr")
   for coord in mse[-1].coords():
     if coord.name() in ['longitude','latitude']:
       coord.long_name=None
 mse = mse.merge_cube()
 if mse.ndim == 2:
   rmse = mse.collapsed("year",iris.analysis.MEAN)**0.5
 else:
   rmse = mse**0.5
 return mse,rmse

"""

def load_all(years12,years2,scratchpath=None,regen=False):
  # load sst data (monthly means)
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
  for cube in sst12+sstref:
    cube.remove_coord("month")
  sst12=sst12.concatenate_cube()
  sst2=sst2.concatenate_cube()
  #sst12=sst12.concatenate_cube().aggregated_by("month",iris.analysis.MIN).collapsed("time",iris.analysis.MEAN)
  sstref=sstref.concatenate_cube()
  sstref.data = np.ma.masked_values(sstref.data,0)
  sst12.data.mask += sst12.data<=2
  sst2.data.mask += sst2.data<=2
  return sst12,sst2,sstref
 

def rmse_new():
  # compute rmse for climatological seasonal mean
  # load data
  MC12 = iris.load_cube("/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/sst/MC12_clim_daily_sea_temperature.nc",cx&cy)
  MC2 =  iris.load_cube("/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/sst/MC2_clim_daily_sea_temperature.nc",cx&cy)
  ref = iris.load_cube("/home/users/emmah/eval/clim_reference_sst.nc",cx&cy)
  ref.coord("longitude").bounds=None
  ref.coord("latitude").bounds=None
  ref=ref.extract(cx&cy)
  ref.remove_coord('time')
  MC2.remove_coord('time')
  ref.add_dim_coord(MC12.coord('time'),0)
  MC2.add_dim_coord(MC12.coord('time'),0)
  # mask (note units are in Kelvin)
  MC12.data.mask += MC12.data < 5
  MC2.data.mask += MC2.data < 5
  ref.data = np.ma.masked_array(ref.data,ref.data==0)
  # compute rmse
  rmse_MC12 = ((MC12 - ref)**2).collapsed(['longitude','latitude'],iris.analysis.MEAN)**0.5
  rmse_MC2 = ((MC2 - ref)**2).collapsed(['longitude','latitude'],iris.analysis.MEAN)**0.5
  return rmse_MC12,rmse_MC2

def main(years12,years2,scratchpath,figname=None,regen=False):
  # create figures
#  mse12, rmse12 = combine_rmse(years12,"MC12",5,scratchpath,regen=regen)
#  mse2, rmse2 = combine_rmse(years2,"MC2",5,scratchpath,regen=regen)
  # load rmse data
  rmse12,rmse2 = rmse_new()
  # load sst data
  sst12,sst2,sstref = load_all(years12,years2,scratchpath,regen)
  if not figname is None:
    ct=iris.Constraint(time=lambda t:t.point.month in [12,1,2])
    # calculate means
    sst12=sst12.collapsed("time",iris.analysis.MEAN)
    sst2=sst2.collapsed("time",iris.analysis.MEAN)
    sstref=sstref.extract(ct).collapsed("time",iris.analysis.MEAN)

    # deal with metadata
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
    # build plots
    fig=bias_plots([sst2[3],sst12[3],sstref[3]], ["MC2","MC12","Reference"] ,cmap1="magma",cmap2="bwr",projection=ccrs.PlateCarree(),minmaxN1rangeN2=(20,32,13,1.1,12),nmax=12,above=True,mark_inner=True)
    fig.set_figwidth(12)
    fig.set_figheight(7)
    fig.subplots_adjust(right=0.94,left=0.04,bottom=0.05)
    fig.suptitle("Foundational Sea Surface Temperature")
    fig.tight_layout()
    # add panel g:rmse
    ax = fig.add_axes([0.07,0.07,0.5,0.2])
    ax.set_xticks([365-31,365-16,365+1,365+15,365+31,365+31+15,365+31+28])
    ax.set_xticklabels(["1-Dec","15-Dec","1-Jan","15-Jan","1-Feb","15-Feb","28-Feb"])
    plt.plot(np.arange(365-31,365-31+90,1),rmse12.data,label="MC12",axes=ax,c="k")
    plt.plot(np.arange(365-31,365-31+90,1),rmse2.data,label="MC2",axes=ax,c="b")
    plt.grid()
    #[iplt.plot(mse12[i]**0.5,axes=ax,label=year,c="0.5",lw=1) for i,year in enumerate(years12)]
    ax.set_title("(g) RMSE Error")
    plt.legend() 
    #ax.set_xlim(365-31,365+31+28)
    fig.savefig(figname)
    plt.show()

if __name__=="__main__":
  main(MC12_years,"/work/scratch-pw2/emmah/eval/","/home/users/emmah/eval_figs/sst_bias.png",regen=False)
