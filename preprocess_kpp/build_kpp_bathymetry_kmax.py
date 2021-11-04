import iris
import numpy as np
from routines import convolve2d
import iris
import matplotlib.pyplot as plt
from scipy.ndimage import label, distance_transform_edt

lsm = iris.load("/gws/nopw/j04/terramaris/emmah/shallow/lsm_ocndepth_mod.nc","lsm")[0]

#bathy = iris.load("/gws/nopw/j04/terramaris/emmah/cmems_ostia/GLO-MFC_001_030_mask_bathy.nc","sea_floor_depth_below_geoid")[0]

cmems = iris.load("/work/scratch-nopw/emmah/cmems/200910/GLOBAL_REANALYSIS_PHY_001_030-TDS_2009-11-01.nc")[0][0]

z = cmems.coord("depth").points

z3d = ((1-cmems.data.mask.T)*z).T

bathy = cmems[0].copy(data=z3d.max(axis=0))
bathy.data=np.ma.masked_array(bathy.data,bathy.data==0)
cx = iris.Constraint(longitude=lambda lon: 85<=lon<=160)
cy = iris.Constraint(latitude=lambda lat: -20<=lat<=20)

bathy.coord("longitude").guess_bounds()
bathy.coord("latitude").guess_bounds()


def interpolate_maximum(x_old,y_old,data,x_new,y_new):
  ix_m = np.array([np.argmin(np.ma.masked_array(np.abs(x_old-x),mask=(x_old>x))) for x in x_new])
  iy_m = np.array([np.argmin(np.ma.masked_array(np.abs(y_old-y),mask=(y_old>y))) for y in y_new])
  ix_p = ix_m+1
  iy_p = iy_m+1
  out = np.ma.masked_array([data[iy_m][:,ix_m],data[iy_m][:,ix_p],data[iy_p][:,ix_m],data[iy_p][:,ix_p]]).max(axis=0)
  return out

plt.subplot(221)
plt.imshow(bathy.extract(cx&cy).data*-1,vmin=-300,vmax=1,origin="lower")

#data=np.maximum(-1*bathy.data,-300)
#data=np.ma.masked_array(data,mask=np.isnan(data))
#bathy=lsm.copy(data = interpolate_maximum(bathy.coord("longitude").points,bathy.coord("latitude").points,data,lsm.coord("longitude").points,lsm.coord("latitude").points))
bathy=bathy.regrid(lsm,iris.analysis.AreaWeighted())

bathy.data =np.maximum(-1*bathy.data,-300)

min_point = -15

z = bathy.data

z=np.ma.masked_array(z,mask=np.isnan(z))
min_coast = z.max()
print(min_coast)

z = np.ma.masked_array(z.filled(0),mask=lsm.data)

z2 = (convolve2d((z.filled(0)>=min_point).astype(float),np.ones((3,3))/9.0)>0.25)+(z.filled(0)>=min_point)


x,n=label(z2)
sizes = np.array([1000]+[(x==i).sum() for i in range(1,n+1)])
a = np.ma.masked_array(sizes[x]<50,lsm.data)
a = a.astype(int)-(x==0).astype(int)
a[z<min_point] = -1
a.mask += lsm.data.astype(bool)


a = (a==1)

x,n = label(a)

znew = z.copy()

ixs,iys = np.meshgrid(range(z.shape[1]),range(z.shape[0]))
x_i,y_i = np.ones(z.shape)*-1,np.ones(z.shape)*-1

plt.subplot(222)
plt.imshow(znew,origin="lower")

for i in range(1,n+1):
  b=distance_transform_edt(x!=i)
  bm = b.min()
  znew[x==i] = z[(b==1)].max()
  y_i[x==i] = iys[b==1][np.argmax(z[b==1])]
  x_i[x==i] = ixs[b==1][np.argmax(z[b==1])]


plt.subplot(223)
plt.imshow(znew,origin="lower")


filled = np.ones(znew.shape)

mask = (z.filled(0)==0)
land_fill  = np.where(znew==0)
for i,j in zip(land_fill[0],land_fill[1]):
  b=np.ma.masked_array(distance_transform_edt(np.hypot(iys-i,ixs-j)),mask=mask)
  bm = b.min()
  znew[i,j] = z[(b==bm)].max()
  y_i[i,j] = iys[b==bm][np.argmax(z[b==bm])]
  x_i[i,j] = ixs[b==bm][np.argmax(z[b==bm])]
  filled[i,j] = z[int(y_i[i,j]),int(x_i[i,j])]

plt.subplot(224)
plt.imshow(znew,origin="lower")

print(znew.max())

plt.show()

plt.imshow(filled,origin="lower")

(znew-z)[znew != z]

out_indices = iris.cube.CubeList([lsm.copy(data=x_i),lsm.copy(y_i)])
out_indices[0].rename("x-indices")
out_indices[1].rename("y-indices")
out = iris.cube.CubeList([lsm,lsm.copy(data=znew.filled(0))])
out[1].rename("max_depth")
out[1].units="m"

iris.save(out,"lsm_ocndepth_kmax.nc")
iris.save(out_indices,"kmax_indices.nc")

