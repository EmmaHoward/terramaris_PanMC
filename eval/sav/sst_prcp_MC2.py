import iris
import iris.cube
import iris.coords
from iris.coord_categorisation import add_day_of_year
import iris.analysis.stats
import numpy as np
dos = iris.coords.DimCoord(np.arange(0,90),units='days')
dos.rename('day')
path="/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/"
tt=iris.cube.CubeList()
pp=iris.cube.CubeList()
for year in [2003,2005]:# [2007,2009,2012,2014,2015,2016,2017,2018]:
    print(year)
    p =iris.load([path+"precip/%04d%02d_hourly_coarsened_rainfall.nc"%(yr,month) for (yr,month) in [(year,12),(year+1,1),(year+1,2)]])
    iris.util.equalise_attributes(p)
    p=p.concatenate_cube()
    t=iris.load_cube(path+"sst/MC2_%d??_daily_sea_temperature.nc"%year)[:,5,:-1,1:-1]
    try:
        add_day_of_year(p,'time','doyr')
    except ValueError:
        pass
    p=p.aggregated_by('doyr',iris.analysis.MEAN)
    p.remove_coord('time')
    p.add_dim_coord(dos,0)
    t.remove_coord('time')
    t.add_dim_coord(dos,0)
    yrc=iris.coords.AuxCoord(year,units='years')
    yrc.rename('Year')
    p.add_aux_coord(yrc)
    t.add_aux_coord(yrc)
    p.remove_coord('doyr')
    t.remove_coord('doyr')
    iris.save(p,"/work/scratch-pw2/emmah/MC2_prcp_%d.nc"%year)
    iris.save(t,"/work/scratch-pw2/emmah/MC2_temp_%d.nc"%year)


"""
    tt.append(t)
    pp.append(p)
    
for cube in tt+pp:
   cube.remove_coord('doyr')

t=tt.merge_cube()
p=pp.merge_cube()

corr = iris.cube.CubeList()
for i in range(20):
    cc = iris.coords.AuxCoord(i,units='days')
    cc.rename('lag')
    cc2 = iris.coords.AuxCoord(-i,units='days')
    cc2.rename('lag')
    if i==0:
        p2 = t[:].copy(data=p[:].data)
        corr.append(iris.analysis.stats.pearsonr(t[:],p2,['day','Year']))
        corr[-1].add_aux_coord(cc)   
    else:     
        p2 = t[:,i:].copy(data=p[:,:-i].data)
        corr.append(iris.analysis.stats.pearsonr(t[:,i:],p2,['day','Year']))
        corr[-1].add_aux_coord(cc)
        p2 = t[:,:-i].copy(data=p[:,i:].data)
        corr.append(iris.analysis.stats.pearsonr(t[:,:-i],p2,['day','Year']))
        corr[-1].add_aux_coord(cc2)

corr=corr.merge_cube()
"""
