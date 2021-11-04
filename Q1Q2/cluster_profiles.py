#!/apps/jasmin/jaspy/miniconda_envs/jaspy3.7/m3-4.9.2/envs/jaspy3.7-m3-4.9.2-r20210320/bin/python
#SBATCH -p test
#SBATCH -n 4
#SBATCH --job-name cluster
#SBATCH -o /home/users/emmah/log/cluster_2.out 
#SBATCH -e /home/users/emmah/log/cluster_2.err 
#SBATCH --time=04:00:00

import datetime as dt
import numpy as np
import iris
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import equalise_attributes
from Q1Q2_routines import read_precip_MC108, read_Q1Q2_MC108, apply_mask 
import os
from sklearn.cluster import KMeans

MC="MC12"
path = {"MC2" :"/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/201516_u-cc339/",
        "MC12":"/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/201516_u-cf309/"}[MC]


def cluster(Q,time,lat,lon,fac,nc,template,z):
  n,nv,nz = Q.shape
  k = KMeans(nc).fit(Q.reshape(n,nv*nz))
  fac_mean = np.ma.masked_array([fac[k.labels_==j].mean(axis=0) for j in range(nc)])
  Q1_c = k.cluster_centers_[:,:nz]/fac_mean
  Q2_c = k.cluster_centers_[:,nz:]/fac_mean
  order = np.argsort((np.maximum(Q1_c[:],0)*z).sum(axis=1)/np.maximum(Q1_c[:],0).sum(axis=1))
  label = k.labels_
  label_cube =np.zeros(template.shape)
  labels=k.labels_
  labels2 = labels.copy()
  Q1_c2 = Q1_c.copy()
  Q2_c2 = Q2_c.copy()
  for j in range(nc):
    b = np.where(order==j)[0][0]
    Q1_c2[j] = Q1_c[b]
    Q2_c2[j] = Q2_c[b]
    label_cube[time[label==b],lat[label==b],lon[label==b]]=j+1
  label_cube = template.copy(data=label_cube)
  label_cube.rename("KMeans label")
  label_cube.units=1
  return label_cube,Q1_c2,Q2_c2


def main(year,nday,threshold,nc,i):
  t0 =dt.datetime(year,12,1)
  dates= [t0+dt.timedelta(i) for i in range(nday)]
  print(dates)
  P = read_precip_MC108(dates,i,path,MC)
  mask = P.data > threshold
  data = read_Q1Q2_MC108(dates,i,path,MC)
#  if not os.path.exists("/work/scratch-nopw/emmah/rho/coarse_rho_%s.nc"%MC):
#    rho = read_rho_MC108(dates,data[0].coord("altitude").points,i,MC)
#    iris.save(rho,"/work/scratch-nopw/emmah/rho/coarse_rho_%s.nc"%MC,zlib=True)
#  else:
#    ct=iris.Constraint(time=lambda t: dt.datetime(t.point.year,t.point.month,t.point.day) in dates)
#    rho=iris.load("/work/scratch-nopw/emmah/rho/coarse_rho_%s.nc"%MC)[0].extract(ct)
#  data.append(rho)
  data.append(P)
  print("all loaded")
  profiles,time,lat,lon= apply_mask(data,mask)
  z= data[0].coord("altitude").points
  dz = np.gradient(z)
  Q = np.array([(profiles["Q1"]*profiles["air_density"]*dz/profiles["P"][:,np.newaxis]).filled(0),(profiles["Q2"]*profiles["air_density"]*dz/profiles["P"][:,np.newaxis]).filled(0)]).transpose([1,0,2])
  fac = (profiles["air_density"]*dz/profiles["P"][:,np.newaxis])
  labels,Q1c,Q2c = cluster(Q,time,lat,lon,fac,nc,P,z)
  Q1c = data.extract("Q1")[0][:nc,:,0,0].copy(data=Q1c)
  Q2c = data.extract("Q2")[0][:nc,:,0,0].copy(data=Q2c)
  for coord in ["longitude","latitude","time"]:
    Q1c.remove_coord(coord)
    Q2c.remove_coord(coord)
  kcoord = iris.coords.DimCoord(range(1,nc+1),units=1)
  kcoord.rename("Kmeans")
  iris.save([Q1c,Q2c,labels],"/work/scratch-nopw/emmah/%s_cluster_labels_%dclusters_%0.02fthreshold.nc"%(MC,nc,threshold))

main(2015,90,0.15,4,9)
