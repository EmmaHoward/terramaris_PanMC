#!/home/users/emmah/miniconda3/envs/iris3/bin/python

from subprocess import call
from time import sleep
import datetime as dt


url_ts = 'python -m motuclient --motu http://my.cmems-du.eu/motu-web/Motu --service-id GLOBAL_REANALYSIS_PHY_001_030-TDS --product-id global-reanalysis-phy-001-030-daily --longitude-min 70 --longitude-max 170 --latitude-min -30 --latitude-max 30 --date-min "{0} 12:00:00" --date-max "{1} 12:00:00" --depth-min 0.493 --depth-max 1062.45 --variable so --variable thetao  --out-dir /work/scratch-nopw/emmah/cmems/{2}/ --out-name GLOBAL_REANALYSIS_PHY_001_030-TDS_{0}.nc     --user ehoward --pwd jerfH23c\$\$'


url_ts = '/home/users/emmah/miniconda3/envs/iris3/bin/python -m motuclient --motu http://my.cmems-du.eu/motu-web/Motu --service-id GLOBAL_REANALYSIS_PHY_001_030-TDS --product-id global-reanalysis-phy-001-030-monthly --longitude-min 70 --longitude-max 170 --latitude-min -30 --latitude-max 30 --date-min "{0} 12:00:00" --date-max "{1} 00:00:00" --depth-min 0.493 --depth-max 1062.45 --variable mlotst  --out-dir /work/scratch-nopw/emmah/cmems/ --out-name GLOBAL_REANALYSIS_PHY_001_030-TDS_{0}.nc     --user ehoward --pwd jerfH23c\$\$'


#yrs = "200304"
n = 5

#call("mkdir -p /work/n02/n02/emmah/cmems/%s/"%yrs,shell=True)
for year in range(2003,2019):
  t1 = dt.datetime(year,12,16)
  t2 = dt.datetime(year+1,2,17)
  call(  url_ts.format("%04d-%02d-%02d"%(t1.year, t1.month,t1.day),"%04d-%02d-%02d"%(t2.year, t2.month,t2.day)),shell=True)


