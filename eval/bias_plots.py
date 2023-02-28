#
# generalised plotting function used in figures 3 and 4
#

import iris
import matplotlib.colors as mcolors
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import iris.plot as iplt
import seaborn as sns
from matplotlib.colors import ListedColormap
from cartopy import feature
def choose_bounds(absmax,absmin,N,steps = [1,2,2.5,5,10]):
    """
    generate round number bounds for colourmaps
    """
    symmetric = False
    if absmax >0 and absmin > 0 and (absmax-absmin)/absmax > 0.9:
      absmin=0
    elif absmax <0 and absmin < 0 and (absmax-absmin)/absmin < -0.9:
      absmax=0    
    elif absmax >0 and absmin < 0:
      symmetric = True
    if symmetric:
      absmax = max(absmax,-absmax)
      absmin = -absmax
    ticks =MaxNLocator(N,symmetric=symmetric,steps=steps)._raw_ticks(absmin,absmax)
    if symmetric:
      delta = ticks[1]-ticks[0]
      ticks = np.concatenate([ticks-delta/2,[ticks[-1]+delta/2]])
    return ticks


def bias_plots(cubes,names,cmap1="viridis",cmap2="bwr",projection=None,minmaxN1rangeN2=None,nmax=12,above=True,mark_inner=True):
  """
  make upper diagonal difference plots between 2 or more datasets
  Parameters
  ==========
     cubes: list of iris cubes to compare
     names: list of labels for each cube. same length as cubes
     cmap1: colormap of diagonal (full field) plots
     cmap2: colormap of off-diagonal (difference field) plots
     projection: cartopy projection for axes
     minmaxN1rangeN2: 5-tuple of vmin, vmax, Nticks for diagonal plots, and then 
                      vmax, Nticks for difference plots (where vmin=-vmax)
                      or None to generate automatic bounds
     nmax: if minmaxN1rangeN2 is none, the maximum number of colorticks
     above: Whether to generate above-diagonal or below-diagonal plots
     mark_inner: whether to mark MC2 domain boundary on all plots
  """
  N = len(cubes)
  assert len(names)==N
  # compute colour bounds
  if minmaxN1rangeN2==None:
    absmax = np.ma.max([cube.collapsed([d.name() for d in cube.dim_coords],iris.analysis.MAX).data for cube in cubes])
    absmin = np.ma.min([cube.collapsed([d.name() for d in cube.dim_coords],iris.analysis.MIN).data for cube in cubes])
    ranges = [[np.ma.abs((cube1-cube2).collapsed([d.name() for d in cube1.dim_coords],iris.analysis.MAX).data) for cube2 in cubes[i+1:]] for (i,cube1) in enumerate(cubes)]  
    absrange = np.ma.max([item for sublist in ranges for item in sublist])
    ticks1 = choose_bounds(absmax,absmin,nmax)
    ticks2 = choose_bounds(absrange,-absrange,nmax)
    N1=len(ticks1)
    N2=len(ticks2)
  else:
    assert len(minmaxN1rangeN2)==5
    vmin,vmax,N1,vrange,N2 = minmaxN1rangeN2
    ticks1 = np.linspace(vmin,vmax,N1)
    ticks2 = np.linspace(-vrange,vrange,N2)
  # generate discrete colormap
  cmap1_ = ListedColormap(sns.color_palette(cmap1,N1-1))
  cmap2_ = ListedColormap(sns.color_palette(cmap2,N2-1))
  # panel counter
  count=0
  fig=plt.figure(figsize=(12,9))
  # iterate through cubes
  for i,cube1 in enumerate(cubes):
    for j,cube2 in enumerate(cubes):
      if i==j:
        # plot diagonal (full fields)
        ax=plt.subplot(N,N,N*i+j+1,projection=projection)  
        ax.set_facecolor("0.5")
        if type(projection) != type(None):
          ax.coastlines()
          ax.set_xlim(cube1.coord("longitude").points.min(),cube1.coord("longitude").points.max())
          ax.set_ylim(cube1.coord("latitude").points.min(),cube1.coord("latitude").points.max())
        iplt.pcolormesh(cube1,cmap=cmap1_,vmin=ticks1[0],vmax=ticks1[-1])
        if mark_inner:
          plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
        if ticks1[0]<0 and ticks1[-1]>0:
          plt.colorbar(ticks=((ticks1[1:]+ticks1[:-1])/2)[::len(ticks1)//5])
        else:
          plt.colorbar(ticks=ticks1[::len(ticks1)//5])
        plt.title("(%s) %s"%(['a','d','f'][i],names[i]))
      if (j>i and above) or (j>i and not above):
        # plot off-diagonal (difference plot)
        ax=plt.subplot(N,N,N*i+j+1,projection=projection)  
        ax.set_facecolor("0.5")
        if type(projection) != type(None):
          ax.coastlines(zorder=11)
          ax.set_xlim(cube1.coord("longitude").points.min(),cube1.coord("longitude").points.max())
          ax.set_ylim(cube1.coord("latitude").points.min(),cube1.coord("latitude").points.max())
        iplt.pcolormesh(cube1-cube2,cmap=cmap2_,vmin=ticks2[0],vmax=ticks2[-1])
        if mark_inner:
          plt.plot([90,155,155,90,90],[-15,-15,15,15,-15],"k",lw=1)
        plt.colorbar()#ticks=((ticks2[1:]+ticks2[:-1])/2)[::len(ticks2)//5])
        plt.title("(%s) %s - %s" %(['b','c','e'][count],names[i],names[j]))
        count +=1
        plt.xlim(85,160)     
        plt.ylim(-20,20)
#        ax.add_feature(feature.LAND,zorder=10,facecolor='0.75')
  return fig

  
  
