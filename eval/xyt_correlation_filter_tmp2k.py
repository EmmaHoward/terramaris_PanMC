from matplotlib.colors import ListedColormap
import seaborn as sns
import iris
from scipy.fftpack import dct
import cftime
from iris.experimental.equalise_cubes import equalise_attributes
#from panMC import panMC
from iris.coord_categorisation import add_day_of_year
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import iris.cube
cy =iris.Constraint(latitude=lambda y: -15<=y<=15)
cx =iris.Constraint(longitude=lambda y: 92<=y<152)
cz = iris.Constraint(pressure=lambda p: p in [200,850])
path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/"
era5path = "/gws/nopw/j04/terramaris/emmah/era5/"
bands = np.array([45,30,18,15,10,6,3])


bins_x = np.array([150,90,50,35,27,22,18,16,14,12,10.5,9.5,8.5,7.5,6.5,5.5,4.5,3.5])
bins_t = np.array([135,55,   35,  25,  20,  16,  14,  12, 10.5, 9.5, 8.5,7.5,6.5,5.5,4.5,3.5,2.5,1.5])
tticks = np.array([30,15,10,5,2])
xticks = np.array([30,20,15,10,5])

def adjust_doyr(cube):
  doyr = cube.coord("doyr").points
  if (cube.coord("year").points%4==0).all():
    doyr[doyr>180] -= 1
  doyr[doyr<180] += 365
  cube.coord("doyr").points = doyr



def band_regress(data1,data2):
  if data1.ndim == 4:
    data1.data=detrend(data1.data.filled(0),axis=1)
    data2.data=detrend(data2.data.filled(0),axis=1)
    F1 = np.fft.fft(dct(dct(data1.data,axis=3,type=2,norm='ortho'),axis=2,type=2,norm='ortho'),axis=1)
    F2 = np.fft.fft(dct(dct(data2.data,axis=3,type=2,norm='ortho'),axis=2,type=2,norm='ortho'),axis=1)
  elif data1.ndim==3:
    data1.data=detrend(data1.data.filled(0),axis=0)
    data2.data=detrend(data2.data.filled(0),axis=0)
    F1 = np.fft.fft(dct(dct([data1.data],axis=3,type=2,norm='ortho'),axis=2,type=2,norm='ortho'),axis=1)
    F2 = np.fft.fft(dct(dct([data2.data],axis=3,type=2,norm='ortho'),axis=2,type=2,norm='ortho'),axis=1)

#  F1x =           dct(dct(data1.data,axis=3,type=2,norm='ortho'),axis=2,type=2,norm='ortho')
#  F2x =           dct(dct(data2.data,axis=3,type=2,norm='ortho'),axis=2,type=2,norm='ortho')
#  F1t = np.fft.fft(data1,axis=1)
#  F2t = np.fft.fft(data2,axis=1)
  nyr,nt,ny,nx = F1.shape
  ft = np.fft.fftfreq(nt,)
  Nx,Ny=240,120
  dx,dy=0.25,0.25
  my = np.arange(0,Ny/2,0.5)/Ny/dy
  mx = np.arange(0,Nx/2,0.5)/Nx/dx
  ftt,myy,mxx = np.meshgrid(np.abs(ft),my,mx,indexing="ij")
  kk = np.hypot(myy,mxx)
  correl = np.zeros((len(bins_t)-1,len(bins_x)-1))
  correlx = np.zeros(len(bins_x)-1)
  correlt = np.zeros(len(bins_t)-1)
  for i in range(len(bins_x)-1):
    print(i)
    for j in range(len(bins_t)-1):
      mask = (1/bins_x[i+1] > kk) * (kk > 1/bins_x[i])* (1/bins_t[j+1] > ftt) * (ftt > 1/bins_t[j])
      G1 = F1*mask
      G2 = F2*mask
      g1 = dct(dct(np.fft.ifft(G1,axis=1),axis=2,type=3,norm='ortho'),axis=3,type=3,norm='ortho')
      g2 = dct(dct(np.fft.ifft(G2,axis=1),axis=2,type=3,norm='ortho'),axis=3,type=3,norm='ortho')
      correl[j,i] = np.corrcoef(g1.real.flatten(),g2.real.flatten())[0,1]
      if i==0:
        mask =  (1/bins_t[j+1] > ftt) * (ftt > 1/bins_t[j])
        G1 = F1*mask
        G2 = F2*mask
        g1 = dct(dct(np.fft.ifft(G1,axis=1),axis=2,type=3,norm='ortho'),axis=3,type=3,norm='ortho')
        g2 = dct(dct(np.fft.ifft(G2,axis=1),axis=2,type=3,norm='ortho'),axis=3,type=3,norm='ortho')
        correlt[j] = np.corrcoef(g1.real.flatten(),g2.real.flatten())[0,1]
    mask =  (1/bins_x[i+1] > kk) * (kk > 1/bins_x[i])
    G1 = F1*mask
    G2 = F2*mask
    g1 = dct(dct(np.fft.ifft(G1,axis=1),axis=2,type=3,norm='ortho'),axis=3,type=3,norm='ortho')
    g2 = dct(dct(np.fft.ifft(G2,axis=1),axis=2,type=3,norm='ortho'),axis=3,type=3,norm='ortho')
    correlx[i] = np.corrcoef(g1.real.flatten(),g2.real.flatten())[0,1]
  return correl,correlx,correlt
     
def read_data(var):
  MC12,obs = iris.cube.CubeList(),iris.cube.CubeList()
  for year in years:
    print(year)
    yc = iris.coords.DimCoord(year,long_name="year",units="years")
    if var=="rainfall":
      data = panMC(year,"MC12","rainfall").load_iris(Constraints=cy,variables=["convective_rainfall_amount","stratiform_rainfall_amount"])
      cube = data[0]+data[1]
      cube = cube.rolling_window("longitude",iris.analysis.MEAN,7)
      ct = iris.Constraint(time = lambda t: (t.point.year,t.point.month) in [(year,12),(year+1,1),(year+1,2)] and not (t.point.month,t.point.day) == (2,29))
      gpm = iris.load(["/gws/nopw/j04/klingaman/datasets/GPM_IMERG/daily/3B-DAY.MS.MRG.3IMERG.%04d.V06.1x1.nc"%year,"/gws/nopw/j04/klingaman/datasets/GPM_IMERG/daily/3B-DAY.MS.MRG.3IMERG.%04d.V06.1x1.nc"%(year+1)],cx&cy&ct).concatenate_cube()
      obs.append(gpm.collapsed("latitude",iris.analysis.MEAN).data)
      lonobs = gpm.coord("longitude").points
      cube = cube.interpolate([("longitude",lonobs)],iris.analysis.Linear())
      MC12.append(cube.data)
      iris.save(cube,scratchpath+"%d_MC12_dailyrain.nc"%year,zlib=True)
      iris.save(obs,scratchpath+"%d_obs_dailyrain.nc"%year,zlib=True)
    elif var=="olr":
        print("implement")
        exit()
    elif var=="eastward_wind":
        ct = iris.Constraint(time=lambda t: t.point.hour in [0,6,12,18])
        cube = panMC(year,"MC12","eastward_wind_pd").load_iris(Constraints=cz&ct)[0]
        add_day_of_year(cube,"time","doyr")
        cube = cube.aggregated_by("doyr",iris.analysis.MEAN)
#
        data = iris.load("/gws/nopw/j04/terramaris/emmah/era5/uv_%04d%02d.nc"%(year,(year+1)%100),cx&cy).extract(var)[0][:90]
        data = data.rolling_window("latitude",iris.analysis.MEAN,2)
        data.add_aux_coord(yc)
        add_day_of_year(data,"time","doyr")
        adjust_doyr(data)
        iris.util.promote_aux_coord_to_dim_coord(data,"doyr")
        data.remove_coord("time")
        data.rename("%s_obs"%var)   
        obs.append(data)
#
        data.coord("longitude").guess_bounds()
        cube.coord("longitude").guess_bounds()
        cube.coord("latitude").guess_bounds()
        cube.coord("longitude").coord_system = data.coord("longitude").coord_system
        cube.coord("latitude").coord_system = data.coord("latitude").coord_system
        cube = cube.regrid(data,iris.analysis.AreaWeighted())
        cube.add_aux_coord(yc)
        adjust_doyr(cube)
        cube.remove_coord("time")
        iris.util.promote_aux_coord_to_dim_coord(cube,"doyr")
        cube.remove_coord("forecast_reference_time")
        cube.rename("%s_MC12"%var)
        MC12.append(cube)
        iris.save([cube,data],"/work/scratch-pw2/emmah/%04d_winds_regridded_MC12.nc"%year)
  iris.util.equalise_attributes(obs)
  MC12 = MC12.merge_cube()
  obs = obs.merge_cube()
  return MC12, obs

#correl_P = main("rainfall")

colours = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple"]
def plot(correls):
  cmap1 = ListedColormap(sns.color_palette("viridis",10))
  fig = plt.figure(figsize=(9,5))
  nx = len(correls)
  ny = 1
  ax1 = fig.add_axes([0.08,0.38,0.12,0.55])
  ax1.grid()
  ax4 = fig.add_axes([0.915,0.38+0.075,0.03,0.4])
  for i,key in enumerate((correls)):
    correl,correlx,correlt = correls[key]
    ax2 = fig.add_axes([0.25+i*0.7/nx,0.1,0.6/nx,0.18])
    ax3 = fig.add_axes([0.25+i*0.7/nx,0.38,0.6/nx,0.55])
    ax2.grid()
    ax3.grid()
    a=ax3.pcolormesh(1/bins_x,1/bins_t,correl,cmap=cmap1,vmin=0,vmax=1)
    ax1.plot(correlt,2/(bins_t[1:]+bins_t[:-1]),c=colours[i])
    ax2.plot(2/(bins_x[1:]+bins_x[:-1]),correlx,c=colours[i])
#  ax1.set_ylim(np.fft.fftshift(freqt)[0],np.fft.fftshift(freqt)[-1])
    ax3.set_yticks(1/tticks)
    ax3.set_yticklabels(["1/%d"%t for t in tticks])
    ax2.set_xticks(1/xticks)
    ax3.set_xticks(1/xticks)
    ax2.set_xticklabels(xticks)
    ax3.set_xticklabels(60//xticks)
    ax2.set_xlim(1/bins_x[0],2/(bins_x[-1]+bins_x[-2]))
    ax3.set_xlim(1/bins_x[0],2/(bins_x[-1]+bins_x[-2]))
    ax3.set_ylim(1/bins_t[0],2/(bins_t[-1]+bins_t[-2]))
    ax2.set_ylim(0,1)
    ax3.set_title("MC12 vs Obs - %s"%key,color=colours[i])
    ax3.set_xlabel("Zonal Wavenumber")
    ax2.set_xlabel("Wavelength (degrees)")
  ax1.set_yticklabels(["1/%d"%t for t in tticks])
  ax1.set_yticks(1/tticks)
  ax1.set_ylim(1/bins_t[0],2/(bins_t[-1]+bins_t[-2])) 
  ax1.set_xlim(0,1)
  ax1.set_ylabel("frequency (1/days)")
  fig.colorbar(a,cax=ax4)

def read_winds(years,var):
  MC,obs = iris.cube.CubeList(),iris.cube.CubeList()
  for year in years:
    yc = iris.coords.DimCoord(year,long_name="year",units="years")
    MC.append( iris.load(path+"/wind/MC2_%04d%02d_winds_regridded0p25.nc"%(year,(year+1)%100),cx&cy).extract("eastward_wind")[0])
    obs.append(iris.load(era5path+"uv_%04d%02d.nc"%(year,(year+1)%100),cx&cy).extract("eastward_wind")[0][:90])
    add_day_of_year(obs[-1],"time","doyr")
    obs[-1].add_aux_coord(yc)
    adjust_doyr(obs[-1])
    obs[-1].remove_coord("time")
    MC[-1].remove_coord("time")
    iris.util.promote_aux_coord_to_dim_coord(obs[-1],"doyr")
    iris.util.promote_aux_coord_to_dim_coord(MC[-1],"doyr")
    try:
      MC[-1].remove_coord("forecast_reference_time")
    except:
      continue
  iris.util.equalise_attributes(obs)
  MC = MC.merge_cube()
  obs = obs.merge_cube()
  obs = obs.rolling_window("latitude",iris.analysis.MEAN,2)
  return MC, obs


scratchpath="/work/scratch-pw2/emmah/eval/"
import pickle
def main(MC2_years,scratchpath,figname):
#  P12,Pobs = read_data("rainfall")
  u2,uobs = read_winds(MC2_years,"eastward_wind")
  import pdb;pdb.set_trace()
#  correl_P,correlx_P,correlt_P = band_regress(P12[:,:90],Pobs[:,:])
#  pickle.dump((correl_P,correlx_P,correlt_P),open(scratchpath+"correl_P.pickle","wb"))
  correl_u850,correlx_u850,correlt_u850 = band_regress(u2[:,1],uobs[:,0])
  pickle.dump((correl_u850,correlx_u850,correlt_u850),open(scratchpath+"MC2_correl_u_850.pickle","wb"))
#  correl_u850,correlx_u850,correlt_u850 = pickle.load(open(scratchpath+"correl_u_850.pickle","rb"))
  correl_u200,correlx_u200,correlt_u200 = band_regress(u2[:90,0],uobs[:,1])
  pickle.dump((correl_u200,correlx_u200,correlt_u200),open(scratchpath+"MC2_correl_u_200.pickle","wb"))
#  correl_u200,correlx_u200,correlt_u200 = pickle.load(open(scratchpath+"correl_u_200.pickle","rb"))
  plot({#"precip":[correl_P,correlx_P,correlt_P],
        "u850":[correl_u850,correlx_u850,correlt_u850],
        "u200":[correl_u200,correlx_u200,correlt_u200] })
  plt.savefig(figname)
  
