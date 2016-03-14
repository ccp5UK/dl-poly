#!/usr/bin/env python3

import statis
import sys

if (len(sys.argv)==2):
  it=int(sys.argv[1])
else:
  print("wrong number of arguments.")    
  sys.exit(2)

print(statis.read('STATIS',it,3))
