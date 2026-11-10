"""Exact tensor Bernstein conversion and subdivision in three variables."""
from fractions import Fraction as F
from math import comb
from itertools import product

def bernstein(poly, box):
 n=tuple(max(e[i] for e in poly) for i in range(3));inds=list(product(*(range(d+1)for d in n)))
 a={i:F(poly.get(i,0)) for i in inds}
 for axis,(left,right) in enumerate(box):
  out={}
  for index in inds:
   j=index[axis];z=list(index)
   value=0
   for q in range(j,n[axis]+1):
    z[axis]=q;value+=a[tuple(z)]*comb(q,j)*left**(q-j)*(right-left)**j
   out[index]=value
  a=out
 for axis in range(3):
  out={}
  for index in inds:
   j=index[axis];z=list(index);value=0
   for q in range(j+1):
    z[axis]=q;value+=a[tuple(z)]*F(comb(j,q),comb(n[axis],q))
   out[index]=value
  a=out
 return n,a

def split(poly,axis):
 n,a=poly
 if n[axis]==0:return poly,poly
 left={};right={};other=[i for i in range(3)if i!=axis]
 for base in product(*(range(n[i]+1)for i in other)):
  z=[0,0,0]
  for i,v in zip(other,base):z[i]=v
  row=[]
  for j in range(n[axis]+1):z[axis]=j;row.append(a[tuple(z)])
  ls=[row[0]];rs=[row[-1]]
  while len(row)>1:
   row=[(x+y)/2 for x,y in zip(row,row[1:])];ls.append(row[0]);rs.append(row[-1])
  for j,(x,y)in enumerate(zip(ls,reversed(rs))):z[axis]=j;left[tuple(z)]=x;right[tuple(z)]=y
 return (n,left),(n,right)

def leaf(bs):
 for i,(n,a) in enumerate(bs):
  m=min(a.values())
  if (m>=0 if i<2 else m>0):return i
 return None

