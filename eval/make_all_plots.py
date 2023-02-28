#
# make all figures for TerraMaris model description paper
#


import matplotlib.pyplot as plt
from subprocess import check_call
import sys
import os
import numpy as np

# functionality to use only select years while model runs were still in progress
years_MC2 = [2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]
years_MC12 = [2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]
years_relax = np.arange(2003,2018)

# path for small intermediate files
scratchpath = "/home/users/emmah/eval/"

# path for figures
figpath = "/home/users/emmah/eval_figs_Feb23/"

check_call("mkdir -p %s"%figpath,shell=True)


# generate figure i
i = int(sys.argv[1])


#for i in [13]:
if 1:
  print(i)
  if i==2: # figure 2: heat salt corrections
    import HSCorr
    HSCorr.plot(scratchpath,figpath+"HSCorr.png")

  if i==3: # figure 3: SST biases
    import sst_bias
    sst_bias.main(years_MC12,years_MC2,scratchpath,figname=figpath+"sst_bias.png")

  if i==4: # removed figure: mixed layer depth bias
    import mld_bias
    mld_bias.main(years_MC12,years_MC2,scratchpath,figname=figpath+"mld_bias.png")

  if i==5: # removed/supplementary figure: surface heat flux bias
    import heatfluxbias
    heatfluxbias.main(years_MC12,years_MC2,figname=figpath+"heatflux_bias.png")

  if i==6: # figure 4: precipitation biases
    import precip_bias
    precip_bias.main(years_MC12,years_MC2,figname=figpath+"precip_bias.png")

  if i==7: # figure 5: diurnal cycle
    import diurnal
    gpm_diurnal.plot(years_MC12,years_MC2,figpath+"diurnal_cycle.png")

  if i==8: # figure 6: mean-state cross sections
    import crossections
    crossections.plot(figpath+"cross_section.png",years_MC12,years_MC2,scratchpath+"/xsection/")

  if i==9: # figure 7:
    import precip_hov
    precip_hov.plot(2015,12,scratchpath,figpath+"precip_201512.png")

  if i==10: # removed figure
    import xyt_correlation_filter
    xyt_correlation_filter.main(years_MC12,years_MC2,scratchpath,figname=figpath+"correl_scales_2.png",calc=False)

  if i==11: # figure 10: maps of standard deviations
     import interannual_std_maps
     interannual_std_maps.main(years_MC12,years_MC2,figpath+"interannual_maps.png")

  if i==12: # figure 11: cross sections of standard deviation
    import interannual_std_vertical
    interannual_std_vertical.plot(years_MC12,years_MC2,scratchpath+"/xsection/",figpath+"interannual_xsection.png")

  if i==13: # figure 12: precip profiles and enso
    fig=plt.figure(figsize=(9,7))
    import precip_interannual
    precip_interannual.plot(years_MC12,years_MC2,figpath+"interannual_precip.png",fig)
    import enso_diff_maps
    enso_diff_maps.plot([2007,2015],figpath+"enso_maps2.png",fig)


  if i==14: # figure 8: propagation of mjo-like and wave-like 
    import lagged_regressions
    lagged_regressions.main(years_MC12,years_MC2,figpath+"propagation_lagregress.png")

  if i==15: # figure 9: mjo rainfall composites
    import mjo_composites
    mjo_composites.plot(scratchpath,figpath+"mjo_composites.png")

  if i==16: # figure 14: mjo propagation in kpp
    import mjo_kpp
    mjo_kpp.main("pcolor",False,scratchpath,years_MC12,years_MC2,figpath+"mjo_kpp_%s.png")

