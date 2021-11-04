import datetime as dt
import numpy as np
import iris
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


def main(year,nday,bins):
  mjo=pd.read_csv('/home/users/emmah/reading/ivar/mjo.csv',index_col=0,names=['x0','x1','mag'],parse_dates=True,dayfirst=True)
  mjo = mjo[dt.datetime(2015,12,1):dt.datetime(2016,3,1)]
  phase = np.ceil(((np.arctan2(mjo.x0,-mjo.x1))/(np.pi/4))%8)*(mjo.mag>=1)
  t0 =dt.datetime(year,12,2)
  dates= [t0+dt.timedelta(i) for i in range(nday)]
  P = read_precip_MC108(dates,9,path,MC)
  cluster = iris.load("/gws/nopw/j04/terramaris/emmah/Q1Q2_analysis/201516_%s_cluster_labels_%dclusters_%0.02fthreshold.nc"%(MC,4,0.15),"kmeans_label")[0][:nday*24]
  data = read_Q1Q2_MC108(dates,9,path,MC)
  if MC=="MC12":
    cube2 = np.zeros(cluster.shape)
    cube2[cluster.data==1] = 3
    cube2[cluster.data==2] = 1
    cube2[cluster.data==3] = 4
    cube2[cluster.data==4] = 2
    cluster.data=cube2
  print("all loaded")
  z= data[0].coord("altitude").points
  data.append(P)
  import pdb;pdb.set_trace()
  for j in range(1,5):
    plt.figure(figsize=(7,7))
    for i,(b0,b1) in enumerate(bins):
      ax1 = plt.subplot(2,len(bins),1+i)
      plt.title("Q1 - [%0.1f,%0.1f]"%(b0,b1))
      a=[-1,1,-1.8,1.8,-1,6.2]
      plt.xlim(a[2*i],a[2*i+1])
      plt.grid()
      ax2 = plt.subplot(2,len(bins),1+i+len(bins))
      a=[-0.5,0.5,-0.85,0.85,-0.2,1.5]
      plt.xlim(a[2*i],a[2*i+1])
      plt.title("Q2 - [%0.1f,%0.1f]"%(b0,b1))
      plt.grid()
      for mjo_phase in range(1,9):
        ct = iris.Constraint(time = lambda t:phase[dt.datetime(t.point.year,t.point.month,t.point.day)]==mjo_phase)
        if cluster.extract(ct) != None:
          mask =(cluster.extract(ct).data==j)
          mask = (P.extract(ct).data>b0)*(P.extract(ct).data<b1)*(cluster.extract(ct).data==j)
          profiles = apply_mask(data.extract(ct),mask,coords = False)
          if len(profiles["Q1"]) > 100:
            ax1.plot(profiles["Q1"].mean(axis=0),z/1000,c[mjo_phase-1])
            ax2.plot(profiles["Q2"].mean(axis=0),z/1000,label=mjo_phase,c=c[mjo_phase-1])
        if i==0 and mjo_phase==8:
          plt.legend()
    plt.savefig("%s_%d.png"%(MC,j))
  plt.show() 


def distribution(year,nday,bins):
  mjo=pd.read_csv('/home/users/emmah/reading/ivar/mjo.csv',index_col=0,names=['x0','x1','mag'],parse_dates=True,dayfirst=True)
  mjo = mjo[dt.datetime(2015,12,1):dt.datetime(2016,2,28)]
  phase = np.ceil(((np.arctan2(mjo.x0,-mjo.x1))/(np.pi/4))%8)*(mjo.mag>=1)
  t0 =dt.datetime(year,12,2)
  dates= [t0+dt.timedelta(i) for i in range(nday)]
#  P = read_precip_MC108(dates,9,path,MC)
#  nt,ny,nx = P.shape
  cluster = iris.load("/work/scratch-nopw/emmah/%s_cluster_labels_%dclusters_%0.02fthreshold.nc"%(MC,4,0.15),"kmeans_label")[0]
  if MC=="MC12":
    cube2 = np.zeros(cluster.shape)
    cube2[cluster.data==1] = 3
    cube2[cluster.data==2] = 1
    cube2[cluster.data==3] = 4
    cube2[cluster.data==4] = 2
    cluster.data=cube2
  nt,ny,nx = cluster.shape
  freq = np.zeros((3,4,8))
  ndays  = [(phase==i).sum() for i in range(1,9)]
  for i,(b0,b1) in enumerate(bins):
    for j in range(1,5):
#      for mjo_phase in range(1,9):
#        ct = iris.Constraint(time = lambda t:phase[dt.datetime(t.point.year,t.point.month,t.point.day)]==mjo_phase)
        if cluster.extract(ct) != None:
          mask = (cluster.extract(ct).data==j)
          #mask = (P.extract(ct).data>b0)*(P.extract(ct).data<b1)*(cluster.extract(ct).data==j)
          freq[i,j-1,mjo_phase-1] = mask.sum()/nx/ny
  freq = freq[:,[0,2,1,3]]
  for i,(b0,b1) in enumerate(bins):
    plt.subplot(1,3,i+1)
    plt.bar(range(1,9),freq[i,0]/ndays/24)
    plt.ylabel("cells/hr")
    plt.title("[%.2f - %.2f]"%(b0,b1))
    plt.xticks([1,2,3,4,5,6,7,8])
    for j in range(1,4):
      plt.bar(range(1,9),freq[i,j]/ndays/24,bottom=np.sum(freq[i,:j]/ndays/24,axis=0))
  plt.show()
  plt.xlabel("MJO Phase")
#  import pdb;pdb.set_trace()
  for j in range(1,5):
    ax=plt.subplot(2,2,j,projection=ccrs.PlateCarree())
    ax.coastlines()
    iplt.pcolormesh(cluster[0].copy(data=(cluster.data==j).mean(axis=0)))
    plt.title(["re-evap","anvil","deep","congestus"][j-1])
    plt.colorbar()


 
bins = [[0.15,1.0],[1.0,5.0],[5.0,20]]
main(2015,35,bins)
#distribution(2015,1,bins)


