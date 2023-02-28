#
# figure 7: hovmoller of rainfall in 2015
# 

import iris
from panMC import panMC
import os
import cmocean
import iris.plot as iplt
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
from routines import cmap_discretise
from iris.coord_categorisation import add_hour,add_day_of_year
from iris.coord_categorisation import add_categorised_coord


def add_day_of_month_float(cube, coord, name='dom'):
    # add iris coordinate for days since 1/month with hours/minutes in decimal
    def _time_of_day(coord, value):
        pt = coord.units.num2date(value)
        return pt.day+pt.hour/24.0+pt.minute/60.0/24.0+pt.second/3600/24.0
    add_categorised_coord(cube, name, coord,
                          _time_of_day)

MC2_path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"

# inner domain constraint
cx = iris.Constraint(longitude=lambda x:90<=x<=155)
cy = iris.Constraint(latitude=lambda y:-15<=y<=15)

names = ["MC2","MC12","GPM"]
bands = ["N","E","S"]

nday={12:31,1:31,2:28}
monname = {12:"December",1:"January",2:"February"}
def calc(year,month,scratchpath):
  # load rainfall for a single month, take latitude-band means, save to precip_hov file
  out = iris.cube.CubeList()
  dates = [dt.datetime(year,month,1)+dt.timedelta(i) for i in range(nday[month])]
  P12 = panMC(year,"MC12","rainfall").load_iris(dates,variables=["stratiform_rainfall_amount","convective_rainfall_amount"])
  P12 = (P12[0]+P12[1]).extract(cx&cy)
  ct = iris.Constraint(time=lambda t: t.point.year==year and t.point.month==month)
  P2=iris.load(MC2_path+"%04d%02d_hourly_coarsened_rainfall.nc"%(year+int(month<6),month),ct)[0]*4
#  P2 = panMC(year,"MC2-tmp","rainfall").load_iris(dates,variables=["stratiform_rainfall_amount"])[0]*4
  TRMM_hr = iris.load(["/gws/nopw/j04/klingaman/datasets/GPM_IMERG/half_hourly/{0}/{0}{1:02d}/3B-HHR.MS.MRG.3IMERG.{0}{1:02d}{2:02d}-*.nc".format(year,month,i) for i in range(1,32)],cx&cy)
  TRMM_hr = TRMM_hr.concatenate_cube()
  for i,(y0,y1) in enumerate([(5,15),(-5,5),(-15,-5)]): 
    cy2 = iris.Constraint(latitude=lambda y: y0<=y<=y1)
    for j,cube in enumerate([P2,P12,TRMM_hr]):
      cube.coord("time").convert_units("days since %04d-%02d-01"%(year,month))
      new = cube.extract(cy2).collapsed("latitude",iris.analysis.MEAN)
      new.rename("%s_%s"%(names[j],bands[i]))
      out.append(new)
  iris.save(out,scratchpath+"precip_hov_%04d%02d.nc"%(year,month))

def plot(year,month,scratchpath,figname):
  # plot hovmollers for year and mont
  if os.path.exists(scratchpath+"precip_hov_%04d%02d.nc"%(year,month)):
      # load cached data
    data=iris.load(scratchpath+"precip_hov_%04d%02d.nc"%(year,month))
  else:
      # calculate cached data
    calc(year,month,scratchpath)
    data=iris.load(scratchpath+"precip_hov_%04d%02d.nc"%(year,month))
  lon_12 = data.extract("MC12_N")[0].coord("longitude").points
  cmap=cmap_discretise(cmocean.cm.rain,12)
  fig,axs=plt.subplots(3,3,figsize=(6.5,8),sharex="all",sharey="all")
  for i,band in enumerate(bands):
    for j,name in enumerate(names):
      # homogenise data
      cube = data.extract("%s_%s"%(name,band))[0]
      if name=="GPM":
        cube = cube.interpolate([("longitude",lon_12)],iris.analysis.Linear())
      if name != "MC2": #MC2 is already aggregated
        add_hour(cube,"time","hour")
        add_day_of_year(cube,"time","doyr")
        cube = cube.aggregated_by(['hour',"doyr"],iris.analysis.MEAN)
      if name == "MC12":
        cube = cube*4 # correct units
      add_day_of_month_float(cube,"time","dom")
      cube.coord("time").convert_units("days since %04d-%02d-01"%(year,month))
      print(cube.coord("time").units)
      # draw plots
      ax=axs[i,j]
      a=iplt.pcolormesh(cube,vmin=0,vmax=6,cmap=cmap,axes=ax,coords=["longitude","dom"])
      ax.set_title("(%s) %s"%('abcdefghi'[i*3++j],name))
      if i==0 and name=='GPM':
         ax.annotate("Melor\n(Nona)",(135,8),xycoords="data",ha='right')
         ax.annotate("TD",(135,15),xycoords="data",ha='left')
      if j==0:
        ax.set_ylabel(["5N - 15N", "5N - 5S", "15S - 5S"][i]+"\n %s %d"%(monname[month],year))
      if i==2:
        ax.set_xlabel("Longitude")
      ax.grid(True)
      ax.set_xticks([100,125,150])
#      ax.set_ylim(0,31)
    fig.subplots_adjust(right=0.8,left=0.15,top=0.95,bottom=0.07)
  cax=fig.add_axes([0.85,0.3,0.05,0.4])
  fig.colorbar(a,cax=cax)
  fig.savefig(figname)

