import iris
import iris.cube

path = "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/postprocessed_outputs/precip/"

path_template = "/gws/nopw/j04/terramaris/emmah/flux_corr_N1280/orog.pp"
template = iris.load(path_template)[0][54:-55,42:-42] # bounds extract inner domain

template.coord('longitude').guess_bounds()
template.coord('latitude').guess_bounds()
 
def regrid(year,month,template):
    data=iris.load(path+"%04d%02d_diurnal_precip.nc"%(year+int(month<6),month))
    regridded = iris.cube.CubeList([])
    for cube in data:
       cube.coord('longitude').guess_bounds()
       cube.coord('latitude').guess_bounds()
       regridded.append(cube.regrid(template,iris.analysis.AreaWeighted()))
    iris.save(regridded,path+"%04d%02d_diurnal_precip_coarsened.nc"%(year+int(month<6),month),zlib=True)

for year in [2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]:
  for month in [12,1,2]:
    print(year,month)
    regrid(year,month,template)


