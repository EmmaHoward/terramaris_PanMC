#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
#SBATCH -p short-serial
#SBATCH -o /home/users/emmah/log/all_%a.o
#SBATCH -e /home/users/emmah/log/all_%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=128000

import panMC
import iris

def vi_kpp_fluxes(MC,year,scratchpath):
  wT_turb = panMC.panMC(year,MC,"wT_turb").load_iris()[0]
  wT_solar = panMC.panMC(year,MC,"wT_solar").load_iris()[0]
  wS_turb = panMC.panMC(year,MC,"wS_turb").load_iris()[0]
  cp = panMC.panMC(year,MC,"cp").load_iris()[0]
  rho = panMC.panMC(year,MC,"sea_water_density").load_iris()[0]
  z = wT_turb.coord('depth')
  int_weights = iris.coords.AuxCoord([0.5]+[1 for i in range(68)]+[0.5],units='1')
  rho_i = (cp*rho).interpolate([('depth',z.points)],iris.analysis.Linear())
  wT_turb_int = iris.analysis.maths.multiply(rho_i*wT_turb*z,int_weights,1).collapsed('depth',iris.analysis.SUM).collapsed('time',iris.analysis.MEAN)
  wT_turb_int.rename(wT_turb.name())
  wT_solar_int = iris.analysis.maths.multiply(rho_i*wT_solar*z,int_weights,1).collapsed('depth',iris.analysis.SUM).collapsed('time',iris.analysis.MEAN)
  wT_solar_int.rename(wT_solar.name())
  wS_turb_int = iris.analysis.maths.multiply(rho_i*wS_turb*z,int_weights,1).collapsed('depth',iris.analysis.SUM).collapsed('time',iris.analysis.MEAN)
  wS_turb_int.rename(wS_turb.name())
  iris.save([wT_turb_int,wT_solar_int,wS_turb_int],scratchpath+"vi_kpp_fluxes_%s_%d.nc"%(MC,year),zlib=True)



vi_kpp_fluxes("MC2",2012,"/home/users/emmah/eval/kpp_fluxes/")
