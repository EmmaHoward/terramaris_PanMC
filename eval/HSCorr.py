#
# plots figure 2: KPP heat and salt balances compared to driving model advection
# outputs of "compute" function are provided in HSCorr.nc. These are the inputs to "plot" function
#

import cartopy.crs as ccrs
import iris
from panMC import panMC
import iris.cube
import surface_heat_flux_bias
import numpy as np
import datetime as dt
import os
from calculus_centred import differentiate as D
from iris.analysis.maths import multiply as mul
import matplotlib.pyplot as plt
import iris.plot as iplt

cmems_path = "/gws/nopw/j04/terramaris/emmah/monthly_cmems/"
era5_path = "/gws/nopw/j04/terramaris/emmah/era5/"
relaxpath = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/"
years = range(2003,2018)

def ideg_to_m(cube):
  # compute inverse of grid spacing in metres
  r = 6371000#cube.coord('latitude').coord_system.semi_major_axis
  y=cube.coord('latitude').points
  idx = iris.coords.AuxCoord(1./(r*np.cos(y*np.pi/180.0)),long_name='x map factor',units='m-1.rad')
  idy = iris.coords.AuxCoord(1./r,long_name='y map factor',units='m-1.rad')
  return idx,idy

# kpp vertical grid spacing
h = iris.coords.AuxCoord([ 1.278019, 1.29535, 1.313157, 1.331461, 1.350283, 1.369644, 1.389568,
    1.410081, 1.431208, 1.452978, 1.475421, 1.498568, 1.522452, 1.547111,
    1.572581, 1.598904, 1.626123, 1.654285, 1.683439, 1.71364, 1.744944,
    1.777413, 1.811113, 1.846116, 1.882499, 1.920344, 1.959742, 2.000791,
    2.043596, 2.088272, 2.134946, 2.183754, 2.234845, 2.288385, 2.344553,
    2.403548, 2.465588, 2.530915, 2.599798, 2.672537, 2.749462, 2.830947,
    2.91741, 3.00932, 3.10721, 3.211683, 3.323424, 3.443223, 3.57198,
    3.710742, 3.860719, 4.023331, 4.200244, 4.393431, 4.605245, 4.838517,
    5.096682, 5.38395, 5.705535, 6.067978, 6.479592, 6.951114, 7.496644,
    8.135098, 8.892416, 9.805218, 10.92685, 12.33824, 14.16832 ],units="m")  # KPP grid  spacing

# kpp vertical coordinate
z = iris.coords.AuxCoord([
    -0.6390094, -1.925694, -3.229947, -4.552257, -5.893128, -7.253092,
    -8.632697, -10.03252, -11.45317, -12.89526, -14.35946, -15.84645,
    -17.35696, -18.89175, -20.45159, -22.03733, -23.64985, -25.29005,
    -26.95892, -28.65745, -30.38675, -32.14793, -33.94219, -35.7708,
    -37.63511, -39.53653, -41.47657, -43.45684, -45.47903, -47.54497,
    -49.65658, -51.81593, -54.02523, -56.28684, -58.60331, -60.97736,
    -63.41193, -65.91018, -68.47553, -71.1117, -73.82271, -76.61291,
    -79.48708, -82.45045, -85.50872, -88.66817, -91.93572, -95.31905,
    -98.82664, -102.468, -106.2537, -110.1958, -114.3075, -118.6044,
    -123.1037, -127.8256, -132.7932, -138.0335, -143.5783, -149.465,
    -155.7388, -162.4542, -169.6781, -177.4939, -186.0077, -195.3565,
    -205.7225, -217.3551, -230.6084],units="m")  # KPP grid

def cmems_mean_adv(year,path):
  # approximate heat and salt advection from nemo reanalysis mean state
  cz = iris.Constraint(depth=lambda d: d<=500)
  # global constants
  rho_w=iris.coords.AuxCoord(997,units="kg/m3")
  cp_w= iris.coords.AuxCoord(4200,units="J/kg/K")
  print(year)
  # load cmems data
  data=iris.load(path+"cmems_mom_uv_%04d.nc"%year,cz)
  new = iris.cube.CubeList()
  u = data.extract("eastward_sea_water_velocity")[0]
  v = data.extract("northward_sea_water_velocity")[0]
  T = data.extract("sea_water_potential_temperature")[0]
  S = data.extract("sea_water_salinity")[0]
  # compute inverse grid spacing
  idx,idy=ideg_to_m(T[:,:,1:-1,1:-1])
  # divergence of currents fields
  div = mul(D(u[:,:,1:-1],"longitude"),idx,2)+D(v[:,:,:,1:-1],"latitude")*idy
  div.rename("divergence")
  zc = u.coord("depth").points
  # compute vertical motion from cotinuity equation
  w = [np.trapz(div.data[:,:i],zc[:i],axis=1) for i in range(1,1+len(zc))]
  w = np.ma.masked_array(w,div.data.mask)
  w = (div*div.coord("depth")).copy(data=w.transpose(1,0,2,3))
  w.rename("vertical_velocity")
  # compute advection terms (u*spatial derivatives of T,S)
  udS  = mul(u[:,1:-1,1:-1,1:-1]*D(S[:,1:-1,1:-1],"longitude"),idx,2)
  vdS  =     v[:,1:-1,1:-1,1:-1]*D(S[:,1:-1,:,1:-1],"latitude")*idy
  wdS  =  -1*w[:,1:-1,:,:]*D(S[:,:,1:-1,1:-1],"depth")
  udT  = mul(u[:,1:-1,1:-1,1:-1]*D(T[:,1:-1,1:-1],"longitude"),idx,2)
  vdT  =     v[:,1:-1,1:-1,1:-1]*D(T[:,1:-1,:,1:-1],"latitude")*idy
  wdT  =  -1*w[:,1:-1,:,:]*D(T[:,:,1:-1,1:-1],"depth")
  Tadv = (udT+vdT+wdT)#.collapsed("time",iris.analysis.MEAN)
  Tadv = Tadv*rho_w*cp_w
  # set metadata
  del Tadv.coord('time').attributes["valid_min"]
  del Tadv.coord('time').attributes["valid_max"]
  Tadv.convert_units("kg.s-3.m-1")
  Sadv = (udS+vdS+wdS)#.collapsed("time",iris.analysis.MEAN)
  Sadv = Sadv*rho_w
  del Sadv.coord('time').attributes["valid_min"]
  del Sadv.coord('time').attributes["valid_max"]
  Sadv.convert_units("kg.m-3.s-1")
  nt,nz,ny,nx =Sadv.shape
  Tadv = Tadv.interpolate([("depth",-z.points)],iris.analysis.Linear(extrapolation_mode="extrapolate"))
  Sadv = Sadv.interpolate([("depth",-z.points)],iris.analysis.Linear(extrapolation_mode="extrapolate"))
  # compute integrals
  Sadv_z = -1*iris.analysis.maths.multiply(Sadv,h,1).collapsed("depth",iris.analysis.SUM).collapsed("time",iris.analysis.MEAN)
  Tadv_z = -1*iris.analysis.maths.multiply(Tadv,h,1).collapsed("depth",iris.analysis.SUM).collapsed("time",iris.analysis.MEAN)
  Sadv_z.rename("cmems salinity advection")
  Tadv_z.rename("cmems heat advection")
  iris.save([Tadv_z,Sadv_z],path+"cmems_mean_adv_%d_neg.nc"%year)
  return Tadv_z,Sadv_z 

def read_hcorr0():
  # read heat corrections from first ocean model level
  hcorr = iris.load(relaxpath+"../HS_corrections/00*_tm_15yr_hcorr.nc")
  hcorr = iris.cube.CubeList([cube[:6] for cube in hcorr])
  iris.util.equalise_attributes(hcorr)
  for cube in hcorr:
    try:
      cube.coord("t").rename("time")
    except:
      1
    cube.coord("latitude").attributes={}
    cube.coord("longitude").attributes={}
    cube.attributes={}
    cube.cell_methods=()
    cube.coord("time").long_name=None
    cube.coord("time").var_name=None
  hcorr = hcorr.concatenate_cube().extract(iris.Constraint(time=lambda t: t.point.month in [12,1,2]))[:,0].collapsed("time",iris.analysis.MEAN)
  return hcorr


def read_corr():
  # read heat and salt corrections
  hcorr = iris.load(relaxpath+"../HS_corrections/00*_tm_15yr_hcorr.nc")
  scorr = iris.load(relaxpath+"../HS_corrections/00*_tm_15yr_scorr.nc")
  hcorr = iris.cube.CubeList([cube[:6] for cube in hcorr])
  scorr = iris.cube.CubeList([cube[:6] for cube in scorr])
  # ocean constants
  rho = panMC(2003,"MC12","sea_water_density").load_iris([dt.datetime(2003,12,1)])[0].collapsed("time",iris.analysis.MEAN)
  cp = panMC(2003,"MC12","cp").load_iris([dt.datetime(2003,12,1)])[0].collapsed("time",iris.analysis.MEAN)
  cp.coords()[0].rename("z")
  rho.coords()[0].rename("z")
  # fix metadata
  iris.util.equalise_attributes(hcorr+scorr)
  for cube in hcorr+scorr:
    try:
      cube.coord("t").rename("time")
    except:
      1
    cube.coord("latitude").attributes={}
    cube.coord("longitude").attributes={}
    cube.attributes={}
    cube.cell_methods=()
    cube.coord("time").long_name=None
    cube.coord("time").var_name=None
  hcorr = hcorr.concatenate_cube().extract(iris.Constraint(time=lambda t: t.point.month in [12,1,2]))
  scorr = scorr.concatenate_cube().extract(iris.Constraint(time=lambda t: t.point.month in [12,1,2]))
  rho.remove_coord("time")
  cp.remove_coord("time")
  nt,nz,ny,nx=hcorr.shape
  # apply dimension factors
  hcorr = hcorr#/rho/cp
  scorr = scorr*rho
  # compute integrals
  hcorr_z = iris.analysis.maths.multiply(hcorr[:,:-1],h,1).collapsed("z",iris.analysis.SUM).collapsed("time",iris.analysis.MEAN)
  scorr_z = iris.analysis.maths.multiply(scorr[:,:-1],h,1).collapsed("z",iris.analysis.SUM).collapsed("time",iris.analysis.MEAN)
  return hcorr_z,scorr_z


def compute(years,scratchpath):
  # calculate data for figure 2 
  hcorr,scorr = read_corr() # load heat and salt corrections
  Q_ref = iris.cube.CubeList()
  Q_relax = iris.cube.CubeList()
  cmems = iris.cube.CubeList()
  # load data
  for year in years:
    Q_ref.append(surface_heat_flux_bias.load_ref(year,era5_path)) # reference heat flux
    Q_relax.append(surface_heat_flux_bias.load_relax(year,relaxpath)) # relaxation run heat flux
    cmems += cmems_mean_adv(year,cmems_path) # cmems T and S advection approximations
    Q_relax[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    Q_ref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
  # fix metadata and compute means
  iris.util.equalise_attributes(Q_ref)
  Q_ref=Q_ref.merge_cube().collapsed("time",iris.analysis.MEAN)
  Q_relax=Q_relax.merge_cube().collapsed("time",iris.analysis.MEAN)
  iris.util.equalise_attributes(cmems)
  cmems = cmems.merge()
  Sadv = cmems.extract("cmems salinity advection")[0].collapsed("time",iris.analysis.MEAN)
  Tadv = cmems.extract("cmems heat advection")[0].collapsed("time",iris.analysis.MEAN)
  mask = np.abs(np.gradient(Tadv.data,axis=0))+np.abs(np.gradient(Tadv.data,axis=1))  >1000
  Tadv.data.mask += mask
  Sadv.data.mask += mask
  for cube in [Sadv, Tadv,hcorr,scorr,Q_ref,Q_relax]:
    cube.coord("latitude").coord_system = None
    cube.coord("longitude").coord_system = None
    try:
      cube.coord("longitude").guess_bounds()
      cube.coord("latitude").guess_bounds()
    except:
      continue
  # regrid to heat and salt correction grid
  Sadv = Sadv.regrid(scorr,iris.analysis.AreaWeighted())
  Tadv = Tadv.regrid(hcorr,iris.analysis.AreaWeighted())
  Q_ref = Q_ref.regrid(hcorr,iris.analysis.AreaWeighted())
  # fix units
  Sadv.convert_units(scorr.units)
  hcorr.convert_units("W/m2")
  hcorr.rename("heat_correction")
  scorr.rename("salinity_correction")
  Q_ref.rename("era5 Q")
  Q_relax.rename("relax Q")
  Tadv.convert_units("W/m2")
  # save to file
  iris.save(iris.cube.CubeList([scorr, hcorr, Q_ref,Q_relax, Sadv, Tadv]),scratchpath+"HScorr.nc")


def plot(scratchpath,figname):
  # load and plot data for figure 2
  # core constants
  rho_w=iris.coords.AuxCoord(997,units="kg/m3")
  cp_w= iris.coords.AuxCoord(4200,units="J/kg/K")
  z_c = iris.coords.AuxCoord(318.1274,units="m")
  z_12 = iris.coords.AuxCoord(237.6925,units="m")
  # load data
  data = iris.load(scratchpath+"/HScorr.nc") 
  Sadv = data.extract("cmems salinity advection")[0]#/z_c/rho_w
  Tadv = data.extract("cmems heat advection")[0]#/z_c/rho_w/cp_w
  hfb = (data.extract("relax Q")[0] - data.extract("era5 Q")[0])
  hcorr = data.extract("heat_correction")[0]#/z_12/rho_w/cp_w
  scorr = data.extract("salinity_correction")[0]#/z_12/rho_w
  # normalise units
  Tadv.convert_units("W/m2")
  hcorr.convert_units("W/m2")
  hfb.convert_units("W/m2")
  scorr.convert_units("g/m2/hour")
  Sadv.convert_units("g/m2/hour")
  fig=plt.figure(figsize=(9,5))
  # plot figure
  for i,cube in enumerate([ scorr, Sadv,  hcorr, Tadv,hfb, hcorr+hfb]):
    ax=plt.subplot(3,2,1+i,projection=ccrs.PlateCarree())
    ax.set_facecolor("0.5")
    ax.coastlines()
    ax.set_title(["(a) Salinity Correction","(b) Salinity Advection", "(c) Heat Correction", "(d) Heat Advection","(e) Heat Flux Bias", "(c)+(e)"][i])
    ax.set_xlim(85,160)
    ax.set_ylim(-20,20)
    ax.plot([90,90,155,155,90],[-15,15,15,-15,-15],"k",lw=1)
    if i in [0,1]:
      iplt.contourf(cube,np.arange(-55,60,10), cmap="PiYG_r",extend="both")
      c=plt.colorbar(ticks=np.arange(-50,60,20))
      c.set_label("g/m2/hr")
    elif i==4:
      iplt.contourf(cube,np.arange(-70,1,10),cmap="Blues_r",extend="max")
      c=plt.colorbar()#ticks=np.arange(-300,301,100))
      c.set_label("W/m2")
    else:
      iplt.contourf(cube,np.arange(-325,326,50),cmap="bwr",extend="both")
      c=plt.colorbar(ticks=np.arange(-300,301,100))
      c.set_label("W/m2")
  fig.tight_layout()
  plt.savefig(figname)
  plt.show()


if __name__=="__main__":
  #compute(range(2003,2018),"/work/scratch-pw2/emmah/eval/")
  plot("/work/scratch-pw2/emmah/eval/","/home/users/emmah/eval_figs/HSCorr.png")
