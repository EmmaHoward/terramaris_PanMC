import iris
from iris.coord_categorisation import add_day_of_year
import iris.analysis.stats
import numpy as np
dos = iris.coords.DimCoord(np.arange(0,90),units='days')
dos.rename('day')
template=iris.load_cube("/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/t_p_corr.nc")
template.coord('latitude').guess_bounds()
template.coord('longitude').guess_bounds()
tt=iris.cube.CubeList()
pp=iris.cube.CubeList()
cx = iris.Constraint(longitude=lambda x:90<=x<=155)
cy = iris.Constraint(latitude=lambda x:-15<=x<=15)
ct = iris.Constraint(time=lambda t: t.point.month in [12,1,2] and not (t.point.month==2 and t.point.day==29))
import os
for year in range(2007,2020):
#  if not os.path.exists("/work/scratch-pw2/emmah/obs_prcp_%d.nc"%year):
    print(year)
    p = iris.load(["/gws/nopw/j04/klingaman/datasets/GPM_IMERG/daily/{year}/3B-DAY.MS.MRG.3IMERG.{year}{month:02}*.nc4".format(year=yr,month=month) for (yr,month) in [(year,12),(year+1,1),(year+1,2)]],cx&cy&ct).extract("Daily accumulated High Quality precipitation from all available MW sources")
#    p=iris.load("/work/scratch-pw2/dship/obs/GPM-IMERG-v6/%d/*.nc"%year,cx&cy&ct)
    iris.util.equalise_attributes(p)
    p=  p.merge_cube()
#    try:
#        add_day_of_year(p,'time','doyr')
#    except ValueError:
#        pass
#    p=p.aggregated_by('doyr',iris.analysis.MEAN)
    p.coord('longitude').coord_system=None
    p.coord('latitude').coord_system=None
    t=iris.load_cube("/gws/nopw/j04/terramaris/emmah/sst_products/ostia/%d*_ostia_sst.nc"%year,cx&cy&ct)
    t.coord('latitude').guess_bounds()
    t.coord('longitude').guess_bounds()    
    p.coord('latitude').guess_bounds()
    p.coord('longitude').guess_bounds()
    p.transpose([0,2,1])
    p=p.regrid(template,iris.analysis.AreaWeighted())
    t=t.regrid(template,iris.analysis.AreaWeighted())
    p.remove_coord('time')
#    p.remove_coord('doyr')
    p.add_dim_coord(dos,0)
    t.remove_coord('time')
    t.add_dim_coord(dos,0)
    yrc=iris.coords.AuxCoord(year,units='years')
    yrc.rename('Year')
    p.add_aux_coord(yrc)
    t.add_aux_coord(yrc)
    p = p[:,:-1,1:-1]
    iris.save(p,"/work/scratch-pw2/emmah/obs_prcp_%d.nc"%year)
    exit()
    iris.save(t,"/work/scratch-pw2/emmah/obs_temp_%d.nc"%year)

"""
corr = iris.cube.CubeList()
for i in range(20):
    print(i)
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
