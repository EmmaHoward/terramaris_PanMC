import iris
import os
import numpy as np
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs

names = ["radar_reflectivity","density","potential_temperature","pressure_rho_grid","cloud_fraction","qcf","qcl","rain","graupel","specific_humidity","upward_air_velocity","eastward_wind","northward_wind"]

stream = "pa"

def nsp(i):
  a=int(np.ceil(np.sqrt(i)))
  b=int(np.ceil(i/a))
  return a,b

data=iris.cube.CubeList()
for name in names:
  tmp = iris.load("/work/scratch-nopw/emmah/postproc_2km_fc/%s/tma_2km_KPPcoupled_%s_%s_20151201.nc"%(stream,stream,name))
  data.append([cube for cube in tmp if cube.name() != "surface_altitude"][0])


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
