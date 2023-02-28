#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
#SBATCH -p test
#SBATCH -o /home/users/emmah/log/eval-%a.o
#SBATCH -e /home/users/emmah/log/eval-%a.e 
#SBATCH --mem=128000
#SBATCH --time=04:00:00
import numpy as np

years_MC12 = [2003,2007,2012,2014,2015,2016,2017,2018]
years_MC2 = [2016]
years_relax = np.arange(2003,2018)

scratchpath = "/home/users/emmah/eval/"
#scratchpath = "/work/scratch-pw2/emmah/eval/"
ppout_12 ="/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/"#monthly_diurnal_sst/"

for i in [7]:
  if i==1:
    import HSCorr
    HSCorr.calc(years_relax,scratchpath)

  if i==16:
    import mjo_kpp
    mjo_kpp.main("pcolor",True,scratchpath)

  if 21<=i<28:
    import mjo_composites
    mjo_composites.calc(i-20,years_MC12,scratchpath)

  if i==28:
    import mjo_composites
    mjo_composites.calc(years_MC12,scratchpath)



