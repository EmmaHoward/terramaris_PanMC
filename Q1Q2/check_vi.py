#!/apps/contrib/jaspy/miniconda_envs/jaspy3.7/m3-4.5.11/envs/jaspy3.7-m3-4.5.11-r20181219/bin/python
#  SBATCH -p test
#SBATCH -p short-serial-4hr
#SBATCH -A short4hr
#SBATCH --array=[3-19]
#SBATCH -o /home/users/emmah/log/budget_N1280/budget_15-%a.o
#SBATCH -e /home/users/emmah/log/budget_N1280/budget_15-%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=48000
#SBATCH --dependency=afterok:6814135


import matplotlib.pyplot as plt
from iris.aux_factory import HybridHeightFactory
from scipy.fftpack import dct,idct
import stratify
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np
import datetime as dt
from iris.aux_factory import HybridHeightFactory
from iris.analysis import calculus_centred,calculus
import os
import sys
from Q1Q2_routines import regrid_z_asl,interp_areaweighted,read_rho_MC108,read_Q1Q2_MC108, read_precip_MC108

#t1 = dt.datetime(int(sys.argv[1][:4]),int(sys.argv[1][4:6]),int(sys.argv[1][6:8]))
t1 = dt.datetime(2015,12,1)

if t1.month >6:
  year = t1.year
else:
  year = t1.year - 1


MC="MC12"
path = {"MC2" :"/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/",
        "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/201516_u-cf309/"}[MC]



def main(date,i):
# calculate Q1 and Q2 for datetime date, coarsened over i*i N1280 grid-cells
  dates = [date]
  data = read_Q1Q2_MC108(dates,i,path,MC)
  rho = read_rho_MC108(dates,data[0].coord("altitude").points,i,MC)
  Q1 = data.extract("Q1")[0]
  Q2 = data.extract("Q2")[0]
  pp = read_precip_MC108(dates,i,path,MC)
  cp_Lv = iris.coords.AuxCoord(1005.7/2.501e6,units="J/K/kg*kg/J")
  Q1=Q1*cp_Lv
  Q1.coord("altitude").var_name=None
  Q2.coord("altitude").var_name=None
  Q1_v = np.trapz((Q1*rho).data.filled(0),x=Q1.coord("altitude").points,axis=1)
  Q2_v = np.trapz((Q2*rho).data.filled(0),x=Q2.coord("altitude").points,axis=1)
  plt.subplot(121)
  plt.xlabel(r"integrated $Q_1*c_p/L_v$ (K/hr)")
  plt.grid()
  plt.ylabel(r"precip (mm/hr)")
  plt.plot(Q1_v.flatten(),pp.data.flatten(),".")
  plt.subplot(122)
  plt.grid()
  plt.xlabel(r"integrated $Q_2$ (1/hr)")
  plt.ylabel(r"precip (mm/hr)")
  plt.plot(Q2_v.flatten(),pp.data.flatten(),".")
  plt.show()


if __name__ == "__main__":
    job = 17# int(os.environ["SLURM_ARRAY_TASK_ID"]) -1 
    main(t1+dt.timedelta(job),9)

