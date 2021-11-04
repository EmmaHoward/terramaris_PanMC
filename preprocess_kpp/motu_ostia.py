from subprocess import call
from time import sleep
import datetime as dt


url = '~/miniconda3/envs/iris3/bin/python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001-TDS --product-id METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2 --longitude-min 80 --longitude-max 160 --latitude-min -25 --latitude-max 25 --date-min "{0}-10-23 12:00:00" --date-max "{1}-03-22 12:00:00" --variable analysed_sst --out-dir /gws/nopw/j04/terramaris/emmah/sst_products/ostia/ --out-name {0}1023_{1}0322_ostia_sst.nc --user ehoward --pwd jerfH23c\$\$'

for year1 in [2016,2017,2018]:
  year2=year1+1
  call(url.format(year1,year2),shell=True)


