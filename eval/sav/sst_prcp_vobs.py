import iris
import iris.cube
import iris.coords
from iris.coord_categorisation import add_day_of_year
import iris.analysis.stats
import numpy as np
from scipy.signal import detrend
dos = iris.coords.DimCoord(np.arange(0,90),units='days')
import matplotlib.pyplot as plt
import iris.plot as iplt

dos.rename('day')
cx = iris.Constraint(longitude=lambda x:90<=x<=155)
cy = iris.Constraint(latitude=lambda x:-15<=x<=2)
path="/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/postprocessed_outputs/"
for model in ['MC12','MC2']:
    tt=iris.cube.CubeList()
    pp=iris.cube.CubeList()
    for year in [2007,2009,2012,2014,2015,2016,2017,2018]:
        print(year)
        p = iris.load_cube("/work/scratch-pw2/emmah/%s_prcp_%d.nc"%(model,year))
        t = iris.load_cube("/work/scratch-pw2/emmah/obs_temp_%d.nc"%(year))
        template = t[:,::5,::3]
        t.data.mask = t.data.mask + t.data < 10
        """
        template.coord('longitude').bounds=None
        template.coord('latitude').bounds=None
        template.coord('longitude').guess_bounds()
        template.coord('latitude').guess_bounds()
        if not t.coord('longitude').has_bounds():
            t.coord('longitude').guess_bounds()
            t.coord('latitude').guess_bounds()
        if not p.coord('longitude').has_bounds():
            p.coord('longitude').guess_bounds()
            p.coord('latitude').guess_bounds()
        p.coord('latitude').coord_system = template.coord('latitude').coord_system
        p.coord('longitude').coord_system = template.coord('longitude').coord_system
        t = t.regrid(template,iris.analysis.AreaWeighted(mdtol=0.5))
        p = p.regrid(template,iris.analysis.AreaWeighted(mdtol=0.5))
        """
        mask_t = t.data.mask
        mask_p = p.data.mask
        p.data=np.ma.masked_array(detrend(p.data,axis=0),mask_p)
        t.data=np.ma.masked_array(detrend(t.data,axis=0),mask_t)
        tt.append(t)
        pp.append(p)
    for cube in pp+tt:
       try:
          cube.remove_coord('doyr')
       except:
           pass
       try:
          cube.remove_coord('forecast_reference_time')
       except:
           pass
    iris.util.equalise_attributes(tt+pp)
    t=tt.merge_cube()
    p=pp.merge_cube()
    print([p,t])
    p.coord('longitude').coord_system=None
    p.coord('latitude').coord_system=None
    p = p.regrid(t,iris.analysis.Nearest())
#    if model in ['obs']:
#      t = t[:,:,:-1,1:-1]
#    if model in ['MC12']:
#      p = p[:,:,:-1,1:-1]
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
    iplt.plot(corr.extract(cy).collapsed(['longitude','latitude'],iris.analysis.MEAN),label=model)
    iris.save(corr,"/work/scratch-pw2/emmah/corr_%s_vobs.nc"%(model))

plt.legend()
plt.show()
