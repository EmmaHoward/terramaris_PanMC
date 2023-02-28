#
# plot sst-precip lagged correlations (figure 13)
# data is generated from sst_prcp_corr
#

import iris
import iris.cube
import iris.coords
from iris.coord_categorisation import add_day_of_year
import iris.analysis.stats
import numpy as np
from scipy.signal import detrend
dos = iris.coords.DimCoord(np.arange(0,90),units='days')
import matplotlib.pyplot as plt
import iris.plot as iplt

dos.rename('day')
cx = iris.Constraint(longitude=lambda x:90<=x<=155)
cy = iris.Constraint(latitude=lambda x:-15<=x<=2)
path="/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/"
corr = {}
colours = ['tab:blue','tab:orange','tab:green','tab:orange','tab:green']
lines = ['-','-','-','--','--']
for i,model in enumerate(['obs','MC12','MC2','MC12_vobs','MC2_vobs'][:3])
    corr[model] = iris.load_cube("/work/scratch-pw2/emmah/corr_%s_0p5.nc"%model,cx&cy)
    iplt.plot(corr[model].extract(cy&cx).collapsed(['longitude','latitude'],iris.analysis.MEAN),label=model,c=colours[i],ls=lines[i])

plt.legend()
plt.grid()
plt.ylabel("Correlation")
plt.xlabel("Lag (days) \n SST leads Precip                Precip leads SST ")
plt.title("SST-Precip Lead-lag Correlation")
plt.subplots_adjust(bottom=0.15)
#plt.savefig("leadlag_1.png")
plt.show()

for i,model in enumerate(corr):
  ax=plt.subplot(4,5,i+1)#projection=ccrs.PlateCarree())
  plt.title(model)
  plt.pcolormesh(corr[model].data.max(axis=0),vmin=0,vmax=0.4)
  plt.colorbar()
  ax=plt.subplot(4,5,i+4)#projection=ccrs.PlateCarree())
  plt.pcolormesh(corr[model].data.min(axis=0),vmin=-0.5,vmax=0)
  plt.colorbar()
  ax=plt.subplot(4,5,i+7)#projection=ccrs.PlateCarree())
  plt.pcolormesh(corr[model].data.argmax(axis=0)-19,vmin=-19,vmax=19,cmap='PiYG')
  plt.colorbar()
  ax=plt.subplot(4,5,i+10)#projection=ccrs.PlateCarree())
  plt.pcolormesh(corr[model].data.argmin(axis=0)-19,vmin=-19,vmax=19,cmap='PiYG')
  plt.colorbar()

 
 
