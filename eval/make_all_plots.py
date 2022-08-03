#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.8/m3-4.9.2/envs/jaspy3.8-m3-4.9.2-r20211105/bin/python
#SBATCH -p short-serial
#SBATCH -o /home/users/emmah/log/all_%a.o
#SBATCH -e /home/users/emmah/log/all_%a.e 
#SBATCH -t 04:00:00
#SBATCH --mem=128000
import matplotlib.pyplot as plt
from subprocess import check_call
import sys
import os
import numpy as np
#i =  int(os.environ["SLURM_ARRAY_TASK_ID"])
#years_MC2 = np.array([2007,2014,2015,2016,2017,2018])
#years_MC12 = np.array([2007,2014,2015,2016,2017,2018])
years_MC2 = [2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]
years_MC12 = [2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]
years_relax = np.arange(2003,2018)

scratchpath = "/home/users/emmah/eval/"
figpath = "/home/users/emmah/eval_figs/"

check_call("mkdir -p %s"%figpath,shell=True)
i = int(sys.argv[1])
if 1:
#for i in [16,4,5,6,7,8,9,10,11,12][:1]:
  if i==2: # fin
    import HSCorr
    HSCorr.plot(scratchpath,figpath+"HSCorr.png")

  if i==3: # fin
    import sst_bias
    sst_bias.main(years_MC12,years_MC2,scratchpath,figname=figpath+"sst_bias.png")

  if i==4: # fin
    import mld_bias
    mld_bias.main(years_MC12,years_MC2,scratchpath,figname=figpath+"mld_bias.png")

  if i==5: # fin
    import heatfluxbias
    heatfluxbias.main(years_MC12,years_MC2,figname=figpath+"heatflux_bias.png")

  if i==6: # fin
    import precip_bias
    precip_bias.main(years_MC12,years_MC2,figname=figpath+"precip_bias.png")

  if i==7: # fin
    import diurnal
    diurnal.plot(years_MC12,years_MC2,figpath+"diurnal_cycle.png")

  if i==8: # fin
    import crossections
    crossections.plot(figpath+"cross_section.png",years_MC12,years_MC2,scratchpath+"/xsection/")

  if i==9: # fin
    import precip_hov
    precip_hov.plot(2015,12,scratchpath,figpath+"precip_201512.png")

  if i==10: # add precip 
    import xyt_correlation_filter
    xyt_correlation_filter.main(years_MC12,years_MC2,scratchpath,figname=figpath+"correl_scales_2.png",calc=False)

  if i==11: # fin
     import interannual_std_maps
     interannual_std_maps.main(years_MC12,years_MC2,figpath+"interannual_maps.png")

  if i==12: # fin
    import interannual_std_vertical
    interannual_std_vertical.plot(years_MC12,years_MC2,scratchpath+"/xsection/",figpath+"interannual_xsection.png")

  if i==13: # add MC2 to second when 2 years available
    fig=plt.figure(figsize=(9,7))
    import precip_interannual
    precip_interannual.plot(years_MC12,years_MC2,figpath+"interannual_precip.png",fig)
    import enso_diff_maps
    enso_diff_maps.plot([2007,2015],figpath+"enso_maps.png",fig)


  if i==14: # fin 
    import lagged_regressions
    lagged_regressions.main(years_MC12,years_MC2,figpath+"propagation_lagregress.png")

  if i==15: # add MC2 when more MJO phases are available for MC2 
    import mjo_composites
    mjo_composites.plot(scratchpath,figpath+"mjo_composites.png")

  if i==16: # add MC2 when more MJO phases are available for MC2
    import mjo_kpp
    mjo_kpp.main("pcolor",False,scratchpath,years_MC12,years_MC2,figpath+"mjo_kpp_%s.png")

