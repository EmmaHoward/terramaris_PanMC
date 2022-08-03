import iris
from iris.coord_categorisation import add_hour,add_month,add_categorised_coord

def add_time_of_day(cube, coord, name='hour'):
    def _time_of_day(coord, value):
        pt = coord.units.num2date(value)
        return pt.hour+pt.minute/60.0+pt.second/3600
    add_categorised_coord(cube, name, coord,
                          _time_of_day)
	
for year in [2003,2005,2007,2009,2012,2014,2015,2016,2017,2018]:
  print(year)
  data = iris.load("/work/scratch-pw2/dship/obs/GPM-IMERG-v6/%d/precipitationCal_half-hourly_pt1deg_MC12domain_*.nc"%year)
  iris.util.equalise_attributes(data)
  data=data.concatenate_cube()
  add_time_of_day(data,'time')
  data=data.aggregated_by('hour',iris.analysis.MEAN)
  iris.save(data,"/work/scratch-pw2/emmah/gpm_%d.nc"%year,zlib=True)
