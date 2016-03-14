#!/usr/bin/env python3 

def lines(i):
  li=i//5
  if i%5>0:
    li+=1
  return li

def read(filename,iter,col):
  with open(filename,'r') as f:
    line=f.readline()
    line=f.readline()
    line=f.readline()
    while (line):
      aux=line.split()
      it=int(aux[0])
      num=int(aux[2])
      ll=lines(num)
      if iter==it:
        cl=lines(col)
        cc=col%5
        if cc==0:
          cc=5    
        for l in range(ll):
          line=f.readline()
          if l+1==cl:
            return line.split()[cc-1]
      else:    
        for l in range(ll):
          line=f.readline()
      line=f.readline()
