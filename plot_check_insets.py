import iris
import os
import numpy as np
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs



region = {"pe":"Java","pf":"Bengkulu"}

name = "cloud_fraction"

ax=plt.subplot(111,projection=ccrs.PlateCarree())
ax.coastlines(color="red")
ax.set_xlim(90,155)
ax.set_ylim(-15,15)


for stream in ["pe","pf"]:
  cube = iris.load("/work/scratch-nopw/emmah/postproc_2km_fc/%s/tma_2km_KPPcoupled_%s_%s_20151201.nc"%(stream,region[stream],name),"bulk_cloud_fraction_in_atmosphere_layer")[0]
  nt,nz,ny,nx = cube.shape
  a=iplt.pcolormesh(cube[nt//2,nz//2],axes=ax,vmin=0,vmax=1)
  
plt.colorbar(a)
plt.show()


names = ['cloud_fraction', 'radar_reflectivity', 'density', 'potential_temperature', 'pressure_rho_grid', 'pressure_theta_grid', 'mixing_ratios', 'specific_humidity', 'winds', 'atmos', 'integral_inst', 'integral_agg', 'rainfall', 'surf_inst', 'surf_agg'][11:]

stream = "pf"


def nsp(i):
  a=int(np.ceil(np.sqrt(i)))
  b=int(np.ceil(i/a))
  return a,b

for name in names:
  if name=="radar_reflectivity" and stream=="pf":
    continue
  data = iris.load("/work/scratch-nopw/emmah/postproc_2km_fc/%s/tma_2km_KPPcoupled_%s_%s_20151201.nc"%(stream,region[stream],name))
  data = iris.cube.CubeList([cube for cube in data if cube.name() != "surface_altitude"])
  mx,my = nsp(len(data))
  fig,ax = plt.subplots(mx,my,subplot_kw = {"projection":ccrs.PlateCarree()})
  if len(data)>1:
    ax=ax.flatten()
  else:
    ax = [ax]
  for i,cube in enumerate(data):
    if len(cube.shape)==4:
      nt,nz,ny,nx = cube.shape
      ax[i].coastlines()
      a=iplt.pcolormesh(cube[nt//2,nz//2],axes=ax[i])
      plt.colorbar(a,ax=ax[i])
      ax[i].set_title(cube.name())
    elif len(cube.shape)==3:
      nt,ny,nx = cube.shape
      ax[i].coastlines()
      a=iplt.pcolormesh(cube[nt//2],axes=ax[i])
      plt.colorbar(a,ax=ax[i])
      ax[i].set_title(cube.name())
  fig.suptitle(name)
  plt.show()  
