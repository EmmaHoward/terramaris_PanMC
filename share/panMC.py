import sys
import datetime as dt
import os
import numpy as np
import warnings

class panMC:
  def __init__(self,year,domain,outname,inset=False):
    assert type(year) == int,"year must be an integer" 
    assert (year <=2023) and (year>2003),"year must be between 2003 and 2023" 
    assert domain in ["MC2-tmp","MC2","MC12"],"domain must be MC2 or MC12" 
    assert inset in [False,"Java","Bengkulu"]
    if domain=="MC12":
      assert inset==False,"no insets in MC12 model" 
    assert os.path.exists(panMC.finalpath(domain,year)),"path does not exist: {0}".format(panMC.finalpath(domain,year))
    if outname=="radsim":
      outname="bt_himawari_8_ahi"
    assert outname in panMC.get_outnames(domain,inset=inset),"panMC filename {0} is invalid. \n"+\
                                     "try panMC.panMC.get_outnames(MC2/MC12,stream,inset=True/False) for valid names" 
    self.year=year
    self.yearstr = "{0}{1}".format(year,(year+1)%100)
    self.domain = domain
    self.outname=outname
    self.inset=inset
    self.stream,self.variables = self.get_outname_info()
    self.filename = self.get_filename()

  reinit_step_pa = {"density":4,"potential_temperature":4,"pressure_rho_grid":4,"pressure_theta_grid":4,\
                  "upward_air_velocity":4,"eastward_wind":4,"northward_wind":4,"specific_humidity":4,\
                  "radar_reflectivity":24,"cloud_fraction":24,"qcf":24,"qcl":24,"rain":24,"graupel":24}


  def get_filename(self):
    # construct filename for with generalised date for file category described by PanMC instance
    filestart = {"MC2":"tma_2km_KPPcoupled_","MC2-tmp":"tma_2km_KPPcoupled_", "MC12":"tma_N1280_KPPcoupled_"}[self.domain]
    if self.stream == "pa" and self.domain[:3] == "MC2" and self.reinit_step_pa[self.outname] < 24:
      datestr = "{year:04d}{month:02d}{day:02d}_{hour:02d}"
    else:
      datestr = "{year:04d}{month:02d}{day:02d}"
    if self.outname =="bt_himawari_8_ahi":
      self.filepath = panMC.finalpath(self.domain,self.year)+self.stream+"/"+filestart+self.outname+"-"+datestr+".nc"
    else:
      self.filepath = panMC.finalpath(self.domain,self.year)+self.stream+"/"+filestart+self.stream+"_"+self.outname+"_"+datestr+".nc"
 
  def finalpath(domain,year):
    # get hardcoded file paths for the various PanMC model domains
    #    domain: model domain name. either "MC2" or "MC12". Temporarily contains "MC2-tmp" to distinguished failed 2015-16 run
    if domain=="MC2-tmp":
      return "/gws/nopw/j04/terramaris/emmah/coupled_2km/production_runs/{0}{1}_u-cc339/".format(year,(year+1)%100)
    elif domain=="MC2":
      return "/gws/nopw/j04/terramaris/panMC_um/MC2_RA2T/{0}{1}_u-cc339/".format(year,(year+1)%100)
    elif domain=="MC12":
      return "/gws/nopw/j04/terramaris/panMC_um/MC12_GA7/{0}{1}_u-cf309/".format(year,(year+1)%100)
    else: 
      assert 0, "invalid domain"

  def get_outnames(domain,stream=None,inset=False):
    # get a list of file category names (model data stream, containing variables) for model data streams
    #    domain: model domain name. either "MC2" or "MC12". Temporarily contains "MC2-tmp" to distinguished failed 2015-16 run
    #    stream: model output stream: pa - pf, or kpp or radsim. 
    #    inset:  full domain or 5-minute limited area inset
    outnames = []
    assert stream is None or stream in ["pa","pb","pc","pd","pe","pf","radsim","kpp"]
    if domain in ["MC2-tmp","MC2"]:
      import file_split_2km
      if stream is None and not inset:
        return [key for key in file_split_2km.file_variables_full.keys()]+[key for key in file_split_2km.file_variables_ocean.keys()]+["bt_himawari_8_ahi"]
      elif stream in ["pa","pb","pc","pd"]:
        return [key for key in file_split_2km.file_variables_full.keys() if file_split_2km.file_variables_full[key][0]==stream]
        assert inset==False
      elif stream == "kpp":
        return [key for key in file_split_2km.file_variables_ocean.keys()]
      elif stream == "radsim":
        return "bt_himawari_8_ahi"
      elif stream in ["pe","pf"] or inset:
        return [key for key in file_split_2km.file_variables_pepf.keys()]
    elif domain=="MC12":
      import file_split_N1280
      if stream is None and not inset:
        return [key for key in file_split_N1280.file_variables.keys()]+[key for key in file_split_N1280.file_variables_ocean.keys()]+["bt_himawari_8_ahi"]
      elif stream in ["pa","pb","pc","pd"]:
        return [key for key in file_split_N1280.file_variables.keys() if file_split_N1280.file_variables[key][0]==stream]
        assert inset==False
      elif stream == "kpp":
        return [key for key in file_split_N1280.file_variables_ocean.keys()]
      elif stream == "radsim":
        return "bt_himawari_8_ahi"
      elif stream in ["pe","pf"] or inset:
         assert 0,"no insets in MC12 run" 


  def get_outname_info(self):
    # get postprocessing metadata (model data stream, containing variables) for file category described by PanMC instance
    if self.domain in ["MC2","MC2-tmp"]:
      import file_split_2km
      if self.inset:
        stream = {"Java":"pe","Bengkulu":"pf"}[inset]
        return stream,file_split_2km.file_variables_pepf[self.outname][1]
      elif self.outname in file_split_2km.file_variables_ocean.keys():
        return "kpp",file_split_2km.file_variables_ocean[self.outname]
      elif self.outname in file_split_2km.file_variables_full.keys():
        tmp = file_split_2km.file_variables_full[self.outname]
        return tmp[0],tmp[2]
      elif self.outname == "bt_himawari_8_ahi":
        return "radsim",["brightness_temperature"]
    if self.domain == "MC12":
      import file_split_N1280
      if self.inset:
         assert 0,"no insets in MC12 run" 
      elif self.outname in file_split_N1280.file_variables_ocean.keys():
        return "kpp",file_split_N1280.file_variables_ocean[self.outname]
      elif self.outname in file_split_N1280.file_variables.keys():
        tmp = file_split_N1280.file_variables[self.outname]
        return tmp[0],tmp[2]
      elif self.outname == "bt_himawari_8_ahi":
        return "radsim",["brightness_temperature"]


  def fill_paths_with_dates(self,dates=None,specify_times=False):
    # Generate list of files containing data valid for a list of dates from file category described by PanMC instance
    # with correct file reinitialiations
    #   dates: list of datetime objects valid during model run. 
    #          default=None triggers 1st Dec - 28th Feb
    #   specify_times: True or False. If true, only files with specified start hour are read.
    #                  this is only intended for files with sub-daily reinitialisation
    if specify_times:
      assert not (dates is None), "datetimes must be supplied if specifying times"
      assert self.domain=="MC2" and stream=="pa" and self.reinit_step_pa[self.outname] < 24, "specify times is for large MC2 pa streams with sub-daily reinitialisation"
    if dates is None:
      dates = [dt.datetime(self.year,12,1) + dt.timedelta(i) for i in range(90)]
    if not specify_times:
      assert np.all([d.hour==0 and d.minute==0 for d in dates])
    if self.stream == "pa" and self.domain[:3] == "MC2" and self.reinit_step_pa[self.outname] < 24:
      # 6-hourly files
      reinit = self.reinit_step_pa[self.outname]/24
      if not specify_times:
        dates = [[d + dt.timedelta(i) for i in np.arange(0,24,reinit)] for d in dates]
        dates = [item for sublist in dates for item in sublist]
    elif self.domain=="MC12" and (self.stream in ["pb","pd"] or (self.stream=="radsim" and self.year==2015)):
      # 6-daily files
      t0 = dt.datetime(self.year,12,1) 
      dates = np.unique([t0 + dt.timedelta( 6.0*((d-t0).days//6) ) for d in dates])
    files = [self.filepath.format(year=date.year,month=date.month,day=date.day,hour=date.hour) for date in dates]
    return files

  def load_iris(self,dates=None,Constraints=None,variables=[],specify_times=False,verbose=True):
    # Read iris CubeList containing data valid for a list of dates from file category described by PanMC instance
    #   dates: list of datetime objects valid during model run.
    #          default=None triggers 1st Dec - 28th Feb
    #   variables: list of variables to extract. empty list => all
    #   Constraints: list of iris constraints to apply. None => unconstrained
    #   specify_times: True or False. If true, only files with specified start hour are read.
    #                  this is only intended for files with sub-daily reinitialisation
    import iris
    from iris.experimental.equalise_cubes import equalise_attributes
    files = self.fill_paths_with_dates(dates,specify_times)
    if dates != None and not specify_times:
      ct = iris.Constraint(time = lambda t: (dt.datetime(t.point.year,t.point.month,t.point.day) in dates 
                      or (dt.datetime(t.point.year,t.point.month,t.point.day)-dt.timedelta(1) in dates and t.point.hour==0 and t.point.minute <2)) 
                      and not dt.datetime(t.point.year,t.point.month,t.point.day,t.point.hour,t.point.minute) == dates[0]  )
      Constraints = ct&Constraints
    if len(variables)>0:
      data = iris.load(files,Constraints).extract(variables)
    else:
      data = iris.load(files,Constraints)
    data=data.concatenate().merge()
    for cube in data:
      # warn the user about hour selection, as this may catch someone out otherwise
      if not cube.coord("time").has_bounds() and not self.stream=="kpp" and verbose:
        t = cube.coord("time").units.num2date(cube.coord("time").points)
        warnings.warn("reading instantaneous data times from %02d:%02d to %02d:%02d UTC (inclusive). Set verbose=False to suppress this message"%(t[0].hour,t[0].minute,t[-1].hour,t[-1].minute))
        break
    return data

  def load_xarray(self,dates=None,specify_times=False,**kwargs):
    # Read xarray dataset containing data valid for a list of dates from file category described by panmc instance
    #   dates: list of datetime objects valid during model run.
    #          default=None triggers 1st Dec - 28th Feb
    #   specify_times: True or False. If true, only files with specified start hour are read.
    #                  this is only intended for files with sub-daily reinitialisation
    import xarray
    files = self.fill_paths_with_dates(dates,specify_times)
    dv = ["latitude_longitude", "forecast_period","forecast_reference_time","forecast_period_0","forecast_reference_time_0","fcst_period_bnds"]
    if self.outname in ["atmos","rainfall","surf_agg","integral_agg"]:
      tmp = []
      for f in files:
        tmp.append(xarray.open_dataset(f,chunks="auto",drop_variables = dv,**kwargs))
      times = [key for key in tmp[0].dims.keys() if "time" in key]
      data = []
      for i,time in enumerate(times):
        a = [DS[[ x for x in DS.data_vars.keys() if(time in DS[x].coords)]] for DS in tmp]
        data.append(xarray.concat(a,dim=time))
      return xarray.merge(data)
    else:
      data=xarray.open_mfdataset(files,drop_variables = dv, **kwargs)
    return data
