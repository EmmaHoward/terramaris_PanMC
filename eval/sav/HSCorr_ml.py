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

cx = iris.Constraint(longitude = lambda x: 84 < x < 161)
cy = iris.Constraint(latitude  =  lambda y: -21 < y < 21)


cmems_path = "/gws/nopw/j04/terramaris/emmah/monthly_cmems/"
era5_path = "/gws/nopw/j04/terramaris/emmah/era5/"
relaxpath = "/gws/nopw/j04/terramaris/emmah/coupled_N1280/relaxation_runs/"
years = range(2003,2018)

def ideg_to_m(cube):
  r = 6371000#cube.coord('latitude').coord_system.semi_major_axis
  y=cube.coord('latitude').points
  idx = iris.coords.AuxCoord(1./(r*np.cos(y*np.pi/180.0)),long_name='x map factor',units='m-1.rad')
  idy = iris.coords.AuxCoord(1./r,long_name='y map factor',units='m-1.rad')
  return idx,idy


def cmems_mean_adv(year,path):
  hmix = iris.load("/work/scratch-pw2/emmah/mld.nc",cx&cy).extract("cmems_hmix")[0]
  #if os.path.exists(path+"cmems_mean_adv_%d_ml.nc"%year):
  #  data=iris.load(path+"cmems_mean_adv_%d_ml.nc"%year)
  #  return(data)
  cz = iris.Constraint(depth=lambda d: d<=500)
  rho_w=iris.coords.AuxCoord(997,units="kg/m3")
  cp_w= iris.coords.AuxCoord(4200,units="J/kg/K")
  print(year)
  data=iris.load(path+"cmems_mom_uv_%04d.nc"%year,cz&cx&cy)
  new = iris.cube.CubeList()
  for cube in data:
    new.append(cube[:,:-1])
    new[-1].data.mask += cube[:,1:].data.mask
  data=new
  u = data.extract("eastward_sea_water_velocity")[0]
  v = data.extract("northward_sea_water_velocity")[0]
  T = data.extract("sea_water_potential_temperature")[0]
  S = data.extract("sea_water_salinity")[0]
  idx,idy=ideg_to_m(T[:,:,1:-1,1:-1])
  div = mul(D(u[:,:,1:-1],"longitude"),idx,2)+D(v[:,:,:,1:-1],"latitude")*idy
  div.rename("divergence")
  z = u.coord("depth").points
  w =  [-np.trapz(div.data[:,:i],-z[:i],axis=1) for i in range(1,1+len(z))]
  w = np.ma.masked_array(w,div.data.mask)
  w = (div*div.coord("depth")).copy(data=w.transpose(1,0,2,3))
  w.rename("vertical_velocity")
  udS  = mul(u[:,1:-1,1:-1,1:-1]*D(S[:,1:-1,1:-1],"longitude"),idx,2)
  vdS  =     v[:,1:-1,1:-1,1:-1]*D(S[:,1:-1,:,1:-1],"latitude")*idy
  wdS  =  -1*w[:,1:-1,:,:]*D(S[:,:,1:-1,1:-1],"depth")
  udT  = mul(u[:,1:-1,1:-1,1:-1]*D(T[:,1:-1,1:-1],"longitude"),idx,2)
  vdT  =     v[:,1:-1,1:-1,1:-1]*D(T[:,1:-1,:,1:-1],"latitude")*idy
  wdT  =  -1*w[:,1:-1,:,:]*D(T[:,:,1:-1,1:-1],"depth")
  Tadv = (udT+vdT+wdT)#.collapsed("time",iris.analysis.MEAN)
  Tadv = Tadv#*rho_w*cp_w
  del Tadv.coord('time').attributes["valid_min"]
  del Tadv.coord('time').attributes["valid_max"]
  Tadv.convert_units("K/day")
  Sadv = (udS+vdS+wdS)#.collapsed("time",iris.analysis.MEAN)
  Sadv = Sadv#*rho_w
  del Sadv.coord('time').attributes["valid_min"]
  del Sadv.coord('time').attributes["valid_max"]
  Sadv.convert_units("0.001/day")
  nt,nz,ny,nx =Sadv.shape
  weights=np.ones(nz)
  z = Sadv.coord("depth")
  mask = np.ma.masked_array(np.ones((nz,ny,nx))*(z.points[:,np.newaxis,np.newaxis] > hmix.data[np.newaxis,1:-1,1:-1])*np.gradient(z.points)[:,np.newaxis,np.newaxis],mask=Tadv[0].data.mask)
  mask = Tadv[0].copy(data=mask)
  mask.units = "m"
  Sadv_z = ((Sadv*mask).collapsed("depth",iris.analysis.SUM)/mask.collapsed("depth",iris.analysis.SUM)).collapsed("time",iris.analysis.MEAN)
  Tadv_z = ((Tadv*mask).collapsed("depth",iris.analysis.SUM)/mask.collapsed("depth",iris.analysis.SUM)).collapsed("time",iris.analysis.MEAN)
  import pdb;pdb.set_trace()
  Sadv_z.rename("cmems salinity advection")
  Tadv_z.rename("cmems heat advection")
  iris.save([Tadv_z,Sadv_z],path+"cmems_mean_adv_%d_ml.nc"%year)
  return Tadv_z,Sadv_z 

def read_hcorr0():
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
  hmix = iris.load("/work/scratch-pw2/emmah/mld.nc").extract("mc12_hmix")[0]
  hcorr = iris.load(relaxpath+"../HS_corrections/00*_tm_15yr_hcorr.nc")
  scorr = iris.load(relaxpath+"../HS_corrections/00*_tm_15yr_scorr.nc")
  hcorr = iris.cube.CubeList([cube[:6] for cube in hcorr])
  scorr = iris.cube.CubeList([cube[:6] for cube in scorr])
  rho = panMC(2003,"MC12","sea_water_density").load_iris([dt.datetime(2003,12,1)])[0].collapsed("time",iris.analysis.MEAN)
  cp = panMC(2003,"MC12","cp").load_iris([dt.datetime(2003,12,1)])[0].collapsed("time",iris.analysis.MEAN)
  cp.coords()[0].rename("z")
  rho.coords()[0].rename("z")
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
  hcorr = (hcorr/rho/cp)[:,:-1]
  scorr = scorr[:,:-1]#*rho
  nt,nz,ny,nx=hcorr.shape
  h = np.array([ 1.278019, 1.29535, 1.313157, 1.331461, 1.350283, 1.369644, 1.389568,
    1.410081, 1.431208, 1.452978, 1.475421, 1.498568, 1.522452, 1.547111,
    1.572581, 1.598904, 1.626123, 1.654285, 1.683439, 1.71364, 1.744944,
    1.777413, 1.811113, 1.846116, 1.882499, 1.920344, 1.959742, 2.000791,
    2.043596, 2.088272, 2.134946, 2.183754, 2.234845, 2.288385, 2.344553,
    2.403548, 2.465588, 2.530915, 2.599798, 2.672537, 2.749462, 2.830947,
    2.91741, 3.00932, 3.10721, 3.211683, 3.323424, 3.443223, 3.57198,
    3.710742, 3.860719, 4.023331, 4.200244, 4.393431, 4.605245, 4.838517,
    5.096682, 5.38395, 5.705535, 6.067978, 6.479592, 6.951114, 7.496644,
    8.135098, 8.892416, 9.805218, 10.92685, 12.33824, 14.16832 ])
  z = np.cumsum(h)
  mask = np.ma.masked_array(np.ones((nz,ny,nx))*(z[:,np.newaxis,np.newaxis] > hmix.data[np.newaxis,:,:])*h[:,np.newaxis,np.newaxis],mask=hcorr[0].data.mask)
  mask = hcorr[0].copy(data=mask)
  mask.units = "m"
  hcorr_z = ((hcorr*mask).collapsed("z",iris.analysis.SUM)/mask.collapsed("z",iris.analysis.SUM)).collapsed("time",iris.analysis.MEAN)
  scorr_z = ((scorr*mask).collapsed("z",iris.analysis.SUM)/mask.collapsed("z",iris.analysis.SUM)).collapsed("time",iris.analysis.MEAN)
  hcorr_z.convert_units("K/day")
  scorr_z.convert_units("0.001/day")
  return hcorr_z,scorr_z


def main():
  hcorr,scorr = read_corr()
  Q_ref = iris.cube.CubeList()
  Q_relax = iris.cube.CubeList()
  cmems = iris.cube.CubeList()
  for year in years:
    Q_ref.append(surface_heat_flux_bias.load_ref(year,era5_path))
    Q_relax.append(surface_heat_flux_bias.load_relax(year,relaxpath))
    cmems += cmems_mean_adv(year,cmems_path)
    Q_relax[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
    Q_ref[-1].coord("time").convert_units("hours since 2003-11-01 00:00:00")
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
#  Sadv = Sadv.regrid(scorr,iris.analysis.AreaWeighted())
#  Tadv = Tadv.regrid(hcorr,iris.analysis.AreaWeighted())
  Q_ref = Q_ref.regrid(hcorr,iris.analysis.AreaWeighted())
  hcorr.rename("heat_correction")
  scorr.rename("salinity_correction")
  Q_ref.rename("era5 Q")
  Q_relax.rename("relax Q")
  iris.save(iris.cube.CubeList([scorr, hcorr, Q_ref,Q_relax, Sadv, Tadv]),"/work/scratch-pw2/emmah/HScorr_3.nc")


def plot():
  data = iris.load("/work/scratch-pw2/emmah/HScorr_3.nc") 
  Sadv = data.extract("cmems salinity advection")[0]#.collapsed("time",iris.analysis.MEAN)
  Tadv = data.extract("cmems heat advection")[0]#.collapsed("time",iris.analysis.MEAN)
  #hfb = data.extract("surface_heat_flux_bias")[0] 
  hcorr = -1*data.extract("heat_correction")[0]
  scorr = -1*data.extract("salinity_correction")[0]
  Sadv.convert_units(scorr.units)
  for i,cube in enumerate([scorr, hcorr,  Sadv, Tadv]):
    ax=plt.subplot(2,2,1+i,projection=ccrs.PlateCarree())
    ax.set_facecolor("0.5")
    ax.coastlines()
    ax.set_title(["Salinity Correction","Heat Correction","Salinity Advection", "Heat Advection"][i])
    ax.set_xlim(85,160)
    ax.set_ylim(-20,20)
    ax.plot([90,90,155,155,90],[-15,15,15,-15,-15],"k",lw=1)
    if i in [0,2]:
      iplt.contourf(cube,np.arange(-0.0055,0.0056,0.001), cmap="PiYG_r",extend="both")
      plt.colorbar(ticks=np.arange(-0.005 ,0.0051,0.001))
    else:
      iplt.contourf(cube,np.arange(-0.0275,0.028,0.005),cmap="bwr_r",extend="both")
      plt.colorbar(ticks=np.arange(-0.025,0.026,0.01))
  plt.show()




if __name__ =="__main__":
  plot()
