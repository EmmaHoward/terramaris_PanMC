import numpy as np
import iris
import scipy.ndimage as nd

year=2018
yrstr = "%04d%02d"%(year,(year+1)%100)
def add_bounds(cube,coord):
    if not cube.coord(coord).has_bounds():
        cube.coord(coord).guess_bounds()
    return(cube)

def interp_atmosgrid(cube,regridder=None):
    cube = add_bounds(cube,'longitude')
    cube = add_bounds(cube,'latitude')
    if regridder is None:
        atmos_cube = iris.load('/gws/nopw/j04/terramaris/emmah/coupled_N1280/sav/lsm_ocndepth_for_kpp.tm.ETOPO2v2c.nc')[0]
#        atmos_cube.coord('grid_longitude').rename("longitude")
#        atmos_cube.coord('grid_latitude').rename("latitude")
        atmos_cube.coord('longitude').coord_system=None
        atmos_cube.coord('latitude').coord_system=None
        atmos_cube = add_bounds(atmos_cube,'longitude')
        atmos_cube = add_bounds(atmos_cube,'latitude')
        regridder = iris.analysis.Linear().regridder(cube,atmos_cube)
    cube_out = regridder(cube)
    return(cube_out,regridder)

var_grid=[('temperature','z','temp','ztemp'),('salinity','z','sal','zsal'),('eastward_sea_water_velocity','depth','u','zvel'),('northward_sea_water_velocity','depth','v','zvel')]

out_cubelist = iris.cube.CubeList()
for invar,inz,outvar,outz in var_grid:
    ct=iris.Constraint(t=lambda tt: (tt.point.day==1) and (tt.point.month==11))
    day = [304,305][(year%4)==0]
    if outvar=="temp":
      cube = iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/%04d/%04d_%05d_temperature_relax.nc"%(year,year,day),ct)[0]
    elif outvar=="sal":
      cube = iris.load("/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/%04d/%04d_%05d_salinity_relax.nc"%(year,year,day),ct)[0]
    else:
      cube = iris.load_cube("/work/scratch-nopw/emmah/cmems/%04d%02d/GLOBAL_REANALYSIS_PHY_001_030-TDS_%d-11-01_uv.nc"%(year,(year+1)%100,year),invar)[0]
    #cube = iris.load_cube("/work/n02/n02/emmah/cmems/201516/GLOBAL_REANALYSIS_PHY_001_030-TDS_2015-11-01.nc",invar)[0]
    print(cube)
    lons = cube.coord('longitude').points
    lats = cube.coord('latitude').points
    cube.remove_coord('longitude')
    cube.remove_coord('latitude')
    lon_coord = iris.coords.DimCoord(lons,standard_name='longitude',units='degrees',var_name='longitude')
    lat_coord = iris.coords.DimCoord(lats,standard_name='latitude',units='degrees',var_name='latitude')
    cube.add_dim_coord(lat_coord,1)
    cube.add_dim_coord(lon_coord,2)
    cube.data.mask = cube.data.mask+np.isnan(cube.data)
    if outvar in ['u','v']:
      new = []
      for i in range(cube.shape[0]):
        inv = cube[i].data
        ind = nd.distance_transform_edt(inv.mask,return_distances=False,return_indices=True)
        new.append(inv[tuple(ind)])
      cube.data=new
    interp_cube,regridder = interp_atmosgrid(cube)
    interp_cube_filled = interp_cube#fill_cube(interp_cube,1e-2,itermax=1e3,verbose=True)

    interp_cube_filled.var_name=outvar
    interp_cube_filled.coord(inz).rename(outz)
    if (interp_cube_filled.coord(outz).points > 0).all():
      interp_cube_filled.coord(outz).points = -1.0*interp_cube_filled.coord(outz).points
#    interp_cube_filled = interp_cube_filled.collapsed('time',iris.analysis.MEAN)
    out_cubelist.append(interp_cube_filled)

iris.save(out_cubelist,"/gws/nopw/j04/terramaris/emmah/%d1101_initcond_N1280_kmax.nc"%year)


