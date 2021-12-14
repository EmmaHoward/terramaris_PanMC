
import numpy as np
years_MC2 = [2016]
years_MC12 = [2016]#,2003,2014,2015,2017]
years_relax = np.arange(2003,2018)

scratchpath = "/home/users/emmah/eval/"
figpath = "/home/users/emmah/eval_figs/"


for i in [7]:
  if i==2:
    import HSCorr
    HSCorr.plot(scratchpath,figpath+"HSCorr.png")

  if i==3:
    import sst_bias
    sst_bias.main(years_MC12,years_MC2,scratchpath,figname=figpath+"sst_bias.png")

  if i==5:
    import precip_bias
    precip_bias.main(years_MC12,years_MC2,figname=figpath+"precip_bias.png")

  if i==6:
    import diurnal
    diurnal.plot(years_MC12,years_MC2,figpath+"diurnal_cycle.png")

  if i==7:
    import crossections
    crossections.plot(figpath+"cross_section.png",years_MC12,years_MC2,scratchpath+"/xsection/")

  if i==8:
    import precip_hov
    precip_hov.plot(2016,12,scratchpath,figpath+"precip_201612.png")

  if i==10:
    import precip_interannual
    precip_interannual.plot(years,figpath+"interannual_precip.png")

  if i==11:
    import interannual_std_vertical
    interannual_std_vertical.plot([2003,2014,2016,2017,2018],scratchpath+"/xsection/",figpath+"interannual_xsection.png")

  if i==12:
   # import interannual_std_maps
   # interannual_std_maps.plot(years,figpath+"interannual_maps.png")
    import enso_diff_maps
    enso_diff_maps.plot([2007,2015],figpath+"interannual_maps.png")

  if i==13:
    import lagged_regressions
    lagged_regressions.main(years_MC2,figpath+"propagation_lagregress.png")

  if i==15:
    import mjo_composites
    mjo_composites.plot(scratchpath,figpath+"mjo_composites.png")

  if i==16:
    import mjo_kpp
    mjo_kpp.main("pcolor",False,scratchpath,figpath+"mjo_kpp_%s.png")

