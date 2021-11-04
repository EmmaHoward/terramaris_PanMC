import datetime as dt
from matplotlib.cm import get_cmap
import numpy as np
import iris
import cmocean
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import equalise_attributes
from Q1Q2_routines import read_precip_MC108, read_Q1Q2_MC108,  apply_mask
from sklearn.cluster import KMeans
import pandas as pd
import cartopy.crs as ccrs
import iris.plot as iplt
c= ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:grey","tab:olive",'tab:cyan']

MC="MC2"
path = {"MC2" :"/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/",
        "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/201516_u-cf309/"}[MC]


names= ["decay","anvil","deep","congestus","dry"][::-1]
def main(year,nday,bins,MC,row):
  path = {"MC2" :"/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/",
        "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/201516_u-cf309/"}[MC]
  mjo=pd.read_csv('/home/users/emmah/reading/ivar/mjo.csv',index_col=0,names=['x0','x1','mag'],parse_dates=True,dayfirst=True)
  mjo = mjo[dt.datetime(2015,12,1):dt.datetime(2016,3,1)]
  phase = np.ceil(((np.arctan2(mjo.x0,-mjo.x1))/(np.pi/4))%8)*(mjo.mag>=1)
  t0 =dt.datetime(year,12,1)
  dates= [t0+dt.timedelta(i) for i in range(nday)]
  P = read_precip_MC108(dates,9,path,MC)
  cluster = iris.load("/gws/nopw/j04/terramaris/emmah/Q1Q2_analysis/201516_%s_cluster_labels_%dclusters_%0.02fthreshold.nc"%(MC,4,0.15),"kmeans_label")[0][:nday*24]
  if MC=="MC12":
    cube2 = np.zeros(cluster.shape)
    cube2[cluster.data==1] = 2
    cube2[cluster.data==2] = 1
    cube2[cluster.data==3] = 4
    cube2[cluster.data==4] = 3
    cluster.data=cube2
  if MC=="MC2":
    cube2 = np.zeros(cluster.shape)
    cube2[cluster.data==1] = 4
    cube2[cluster.data==2] = 3
    cube2[cluster.data==3] = 2
    cube2[cluster.data==4] = 1
    cluster.data=cube2
  cluster=cluster.data
  print("all loaded")
  for lag in range(1,7):
    plt.subplot(2,6,lag+6*row,aspect=1)
    t = np.zeros((5,5))
    for i in range(5):
      for j in range(5):
         t[i,j] = ((cluster[lag:]==i)*(cluster[:-lag]==j)).sum()
    if lag == 1:
      plt.yticks([0.5,1.5,2.5,3.5,4.5],names[::-1])
    if row == 1:
      plt.xticks([0.5,1.5,2.5,3.5,4.5],names,rotation=90)
    t = t/(t.sum(axis=1))
    plt.pcolormesh(t[::-1],cmap=cmocean.cm.rain,vmin=0,vmax=1)
    import pdb;pdb.set_trace()

bins = [[0.15,20]]
main(2015,90,bins,"MC2",0)
main(2015,90,bins,"MC12",1)
plt.show()
#distribution(2015,1,bins)


