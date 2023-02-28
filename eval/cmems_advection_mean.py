#
# approximate heat and salt advection from nemo reanalysis mean state
# feeds into figure 2
# 


#!/home/users/emmah/miniconda3/envs/iris3/bin/python
#SBATCH -p short-serial-4hr
#SBATCH -o /home/users/emmah/log/precip.o
#SBATCH -e /home/users/emmah/log/precip.e 
#SBATCH -t 01:00:00
#SBATCH --array=[1-5]

import os
import iris
import numpy as np
from iris.analysis.calculus_centred import differentiate as D
from iris.analysis.maths import multiply as mul
import matplotlib.pyplot as plt
import iris.plot as iplt

def ideg_to_m(cube):
  # compute grid spacing in metres  
  r = 6371000#cube.coord('latitude').coord_system.semi_major_axis
  y=cube.coord('latitude').points
  idx = iris.coords.AuxCoord(1./(r*np.cos(y*np.pi/180.0)),long_name='x map factor',units='m-1.rad')
  idy = iris.coords.AuxCoord(1./r,long_name='y map factor',units='m-1.rad')
  return idx,idy

cz = iris.Constraint(depth=lambda d: d<=350)

rho_w=iris.coords.AuxCoord(997,units="kg/m3")
cp_w= iris.coords.AuxCoord(4200,units="J/kg/K")

#job =  int(os.environ["SLURM_ARRAY_TASK_ID"])-1
#year = [2003,2004,2005,2016,2017][job]

for year in range(2003,2018):
  print(year)
  # load data
  data=iris.load("/gws/nopw/j04/terramaris/emmah/monthly_cmems/cmems_mom_uv_%04d.nc"%year,cz)
  u = data.extract("eastward_sea_water_velocity")[0]
  v = data.extract("northward_sea_water_velocity")[0]
  T = data.extract("sea_water_potential_temperature")[0]
  S = data.extract("sea_water_salinity")[0]
  # grid spacings
  idx,idy=ideg_to_m(T[:,:,1:-1,1:-1])
  # compute current divergence
  div = mul(D(u[:,:,1:-1],"longitude"),idx,2)+D(v[:,:,:,1:-1],"latitude")*idy
  div.rename("divergence")
  z = u.coord("depth").points
  # derive vertical motion from continuity equation
  w =  [-np.trapz(div.data[:,:i],-z[:i],axis=1) for i in range(1,1+len(z))]
  w = (div*div.coord("depth")).copy(data=np.array(w).transpose(1,0,2,3))
  w.rename("vertical_velocity")
  # compute advection terms (u*derivative of S, T)
  udS  = mul(u[:,1:-1,1:-1,1:-1]*D(S[:,1:-1,1:-1],"longitude"),idx,2)
  vdS  =     v[:,1:-1,1:-1,1:-1]*D(S[:,1:-1,:,1:-1],"latitude")*idy
  wdS  =  -1*w[:,1:-1,:,:]*D(S[:,:,1:-1,1:-1],"depth")
  udT  = mul(u[:,1:-1,1:-1,1:-1]*D(T[:,1:-1,1:-1],"longitude"),idx,2)
  vdT  =     v[:,1:-1,1:-1,1:-1]*D(T[:,1:-1,:,1:-1],"latitude")*idy
  wdT  =  -1*w[:,1:-1,:,:]*D(T[:,:,1:-1,1:-1],"depth")
  Tvadv = wdT#.collapsed("time",iris.analysis.MEAN)
  Tadv = (udT+vdT+wdT)#.collapsed("time",iris.analysis.MEAN)
#  Tvadv=((udT[:,0]*u.coord('depth')[0]).copy(data=np.trapz((wdT).data,x=z[1:-1],axis=1)))
#  Tadv=((udT[:,0]*u.coord('depth')[0]).copy(data=np.trapz((udT+vdT+wdT).data,x=z[1:-1],axis=1)))
  Tadv = Tadv*rho_w*cp_w
  Tvadv = Tvadv*rho_w*cp_w
#  adv.append(mul(D(ut_z[:,1:-1],'longitude'),idx[1:-1],1)+ D(vt_z[:,:,1:-1],'latitude')*idy)
  # fix metadata
  del Tvadv.coord('time').attributes["valid_min"]
  del Tvadv.coord('time').attributes["valid_max"]
  del Tadv.coord('time').attributes["valid_min"]
  del Tadv.coord('time').attributes["valid_max"]
  Tadv.rename("total heat advection")
  Tvadv.rename("vertical heat advection")
  Tadv.convert_units("kg.s-3.m-1")
  Tvadv.convert_units("kg.s-3.m-1")
  Svadv = wdS#.collapsed("time",iris.analysis.MEAN)
  Sadv = (udS+vdS+wdS)#.collapsed("time",iris.analysis.MEAN)
#  Svadv=((udS[:,0]*u.coord('depth')[0]).copy(data=np.trapz((wdS).data,x=z[1:-1],axis=1)))
#  Sadv=((udS[:,0]*u.coord('depth')[0]).copy(data=np.trapz((udS+vdS+wdS).data,x=z[1:-1],axis=1)))
  Sadv = Sadv*rho_w
  Svadv = Svadv*rho_w
#  adv.append(mul(D(ut_z[:,1:-1],'longitude'),idx[1:-1],1)+ D(vt_z[:,:,1:-1],'latitude')*idy)
  del Svadv.coord('time').attributes["valid_min"]
  del Svadv.coord('time').attributes["valid_max"]
  del Sadv.coord('time').attributes["valid_min"]
  del Sadv.coord('time').attributes["valid_max"]
  Sadv.rename("total salinity advection")
  Svadv.rename("vertical salinity advection")
  Sadv.convert_units("kg.m-3.s-1")
  Svadv.convert_units("kg.m-3.s-1")
  # save data to file
  iris.save([Tadv,Tvadv,Sadv,Svadv],"/work/scratch-nopw/emmah/coupled/cmems_mean_adv_z_%d.nc"%year)



