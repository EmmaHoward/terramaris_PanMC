import datetime as dt
from matplotlib.cm import get_cmap
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

MC="MC12"
path = {"MC2" :"/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/",
        "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/201516_u-cf309/"}[MC]


names= ["decay","anvil","deep","congestus"]
def main(year,nday,bins):
  mjo=pd.read_csv('/home/users/emmah/reading/ivar/mjo.csv',index_col=0,names=['x0','x1','mag'],parse_dates=True,dayfirst=True)
  mjo = mjo[dt.datetime(2015,12,1):dt.datetime(2016,3,1)]
  phase = np.ceil(((np.arctan2(mjo.x0,-mjo.x1))/(np.pi/4))%8)*(mjo.mag>=1)
  t0 =dt.datetime(year,12,1)
  dates= [t0+dt.timedelta(i) for i in range(nday)]
  P = read_precip_MC108(dates,9,path,MC)
  cluster = iris.load("/gws/nopw/j04/terramaris/emmah/Q1Q2_analysis/201516_%s_cluster_labels_%dclusters_%0.02fthreshold.nc"%(MC,4,0.15),"kmeans_label")[0][:nday*24]
  data = read_Q1Q2_MC108(dates,9,path,MC)
  if MC=="MC12":
    cube2 = np.zeros(cluster.shape)
    cube2[cluster.data==1] = 3
    cube2[cluster.data==2] = 4
    cube2[cluster.data==3] = 1
    cube2[cluster.data==4] = 2
    cluster.data=cube2
  print("all loaded")
  z= data[0].coord("altitude").points
  data.append(P)
  print(data)
  plt.figure(figsize=(7,7))
  colours=  {}
  for j in range(1,5):
      ax1 = plt.subplot(2,4,5-j)
      plt.title("Q1 - %s"%names[j-1])
      plt.grid()
      ax2 = plt.subplot(2,4,9-j)
      plt.title("Q2 - %s"%names[j-1])
      plt.grid()
      for mjo_phase in [8]:
        ct = iris.Constraint(time = lambda t:phase[dt.datetime(t.point.year,t.point.month,t.point.day)]==mjo_phase)
        if cluster.extract(ct) != None:
          mask =(cluster.extract(ct).data==j)*(P.extract(ct).data<1.5)*(P.extract(ct).data>0.5)
          profiles = apply_mask(data.extract(ct),mask,coords = False)
          if len(profiles["Q1"]) > 100:
            rr = np.random.randint(0,len(profiles["Q1"]),30)
            colours[j] = get_cmap("viridis")(profiles["P"][rr]-0.5).copy()
            for i,r in enumerate(rr):
              ax1.plot(profiles["Q1"][r],z/1000,lw=0.5,c=colours[j][i])
              ax2.plot(profiles["Q2"][r],z/1000,lw=0.5,c=colours[j][i])
            ax1.plot(profiles["Q1"].mean(axis=0),z/1000,c="k")#c[mjo_phase-1])
            ax2.plot(profiles["Q2"].mean(axis=0),z/1000,c="k")#label=mjo_phase,c=c[mjo_phase-1])
            plt.draw()
#        if i==0 and mjo_phase==8:
#          plt.legend()
#    plt.savefig("%s_%d.png"%(MC,j))
  plt.show() 
  import pdb;pdb.set_trace()
 
bins = [[0.15,20]]
main(2015,90,bins)
#distribution(2015,1,bins)


