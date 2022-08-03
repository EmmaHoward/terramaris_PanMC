#!/usr/bin/env python

import os
import sys
import numpy as np

prefix = ["B","KiB","MiB","GiB","TiB"]

def getsize(base,name=None):
  out,count = 0,0
  try:
    for key in os.listdir(base):
      filepath = os.path.join(base,key)
      if os.path.isfile(filepath) and not os.path.islink(filepath):
        if name == None or name in filepath:
          out += os.path.getsize(filepath)
          count += 1
      elif os.path.isdir(filepath) and not os.path.islink(filepath):
        tmp = getsize(filepath)
        out += tmp[0]
        count += tmp[1]
  except OSError:
    return out,count
  return out,count

def number_format(x,N):
  if N:
    n = np.where(np.array(prefix)==N)[0][0]
  else: 
    n = int(np.floor(np.log2(x)/10))
    if n>4:
      n=4
  h = prefix[n]
  return x/(1024**n),h

def main(base='.',N=None,name=None):
  print(N)
  if N:
    assert(N in prefix)
  keys = os.listdir(base)
  out = {}
  count = {}
  for key in keys:
    out[key] = 0
    count[key] = 0
    filepath = os.path.join(base,key)
    if os.path.isfile(filepath) and not os.path.islink(filepath):
      if name == None or name in filepath:
        out[key] += os.path.getsize(filepath)
        count[key] += 1
    elif os.path.isdir(filepath) and not os.path.islink(filepath):
      tmp = getsize(filepath,name=name)
      out[key] += tmp[0]
      count[key] += tmp[1]
    if out[key]>0:
      x,h = number_format(out[key],N)
      print("%.2f \t n=%d \t"%(x,count[key]),h,key)
      #print("%.2f \t n=%d \t"%(24*x/count[key],count[key]),h,key)

if len(sys.argv)==4:
  main(sys.argv[1],sys.argv[2],sys.argv[3])
  # Base dir, Prefix, grep
elif len(sys.argv)==3:
  main(sys.argv[1],sys.argv[2])
elif len(sys.argv)==2:
  main(sys.argv[1])
else:
  main('.')
