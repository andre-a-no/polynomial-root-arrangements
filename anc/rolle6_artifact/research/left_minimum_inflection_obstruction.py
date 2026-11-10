"""Exact coupled left-minimum/second-inflection obstruction for real sextics.

The checker uses only integer and rational arithmetic and the standard library.
"""
import json,math,itertools
from fractions import Fraction as F
from pathlib import Path
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]

def add(*ps):
 out={}
 for p in ps:
  for m,c in p.items():out[m]=out.get(m,0)+c
 return {m:c for m,c in out.items()if c}
def scale(p,c):return {m:c*v for m,v in p.items()if c*v}
def mul(p,q):
 out={}
 for m,a in p.items():
  for n,b in q.items():
   r=tuple(x+y for x,y in zip(m,n));out[r]=out.get(r,0)+a*b
 return {m:c for m,c in out.items()if c}
def power(p,n):
 out={(0,0,0):1}
 for _ in range(n):out=mul(out,p)
 return out
def decode(data):
 assert all(len(m)==3 and all(type(i)is int and i>=0 for i in m)and type(c)is int for m,c in data)
 assert len({tuple(m)for m,c in data})==len(data)
 return {tuple(m):c for m,c in data if c}

def raw_inequalities():
 one={(0,0,0):1};K={(1,0,0):1};S={(0,1,0):1};T={(0,0,1):1}
 D=add(power(K,2),scale(K,-1),one)
 AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
 CN=add(scale(power(K,2),2),scale(K,-2),scale(one,-1))
 DN=scale(add(AN,CN),-1);BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
 assert add(mul(AN,CN),mul(AN,DN),mul(CN,DN),scale(power(D,2),3))=={}
 assert add(mul(mul(AN,CN),DN),BN)=={}
 W=add(AN,mul(S,add(scale(D,-1),scale(AN,-1))))
 Z=add(scale(D,-4),mul(T,add(AN,scale(D,4))))
 A4=add(scale(power(W,4),-1),scale(mul(power(D,2),power(W,2)),6),scale(mul(BN,W),-4))
 B5=add(scale(power(Z,5),-1),scale(mul(power(D,2),power(Z,3)),10),scale(mul(BN,power(Z,2)),-10),scale(mul(A4,Z),-5))
 def P(Y):return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20),scale(mul(A4,power(Y,2)),15),scale(mul(B5,Y),6))
 def Q(Y):return add(power(Y,5),scale(mul(power(D,2),power(Y,3)),-10),scale(mul(BN,power(Y,2)),10),scale(mul(A4,Y),5),B5)
 def R(Y):return add(power(Y,4),scale(mul(power(D,2),power(Y,2)),-6),scale(mul(BN,Y),4),A4)
 ys=[AN,scale(D,-1),CN,{},D,DN,W,Z]
 assert R(W)=={} and Q(Z)=={}
 out={}
 for j,sg in enumerate([-1,-1,1,1,-1,-1,-1]):out[f'Q{j}']=scale(Q(ys[j]),sg)
 for j,sg in enumerate([-1,1,1,-1,-1,-1]):out[f'R{j}']=scale(R(ys[j]),sg)
 out['Rleft']=R(Z)
 for i in [0,6]:
  for j in [1,2,3,4,5,7]:out[f'H{i}{j}']=add(P(ys[i]),scale(P(ys[j]),-1))
 out['physical']=scale(CN,-1)
 return out

def bernstein(poly,degrees=None):
 unit={}
 for (i,j,l),c in poly.items():
  for h in range(i+1):unit[h,j,l]=unit.get((h,j,l),F(0))+c*math.comb(i,h)*F(3,8)**h
 unit={m:c for m,c in unit.items()if c}
 natural=tuple(max((m[j]for m in unit),default=0)for j in range(3))
 degrees=natural if degrees is None else tuple(degrees)
 assert all(a<=b for a,b in zip(natural,degrees))
 shape=tuple(n+1 for n in degrees);out=[]
 # Positive change of coordinates from power coefficients to Bernstein.
 for index in itertools.product(*(range(n)for n in shape)):
  value=F(0)
  for m,c in unit.items():
   if all(h<=i for h,i in zip(m,index)):
    weight=F(1)
    for h,i,n in zip(m,index,degrees):weight*=F(math.comb(i,h),math.comb(n,h))
    value+=c*weight
  out.append(value)
 return shape,out

def split(a,shape,axis):
 left=[None]*len(a);right=[None]*len(a)
 strides=(shape[1]*shape[2],shape[2],1);n=shape[axis]-1
 others=[j for j in range(3)if j!=axis]
 for ix in itertools.product(*(range(shape[j])for j in others)):
  base=sum(i*strides[j]for i,j in zip(ix,others));indices=[base+i*strides[axis]for i in range(n+1)]
  row=[a[i]for i in indices];l=[row[0]];r=[row[-1]]
  for _ in range(n):row=[(x+y)/2 for x,y in zip(row,row[1:])];l.append(row[0]);r.append(row[-1])
  for i,x,y in zip(indices,l,reversed(r)):left[i]=x;right[i]=y
 assert all(x is not None for x in left+right)
 return left,right

MASK={0:1,1:-1,4:-1,9:-1,10:1,11:-1}
def reasons(word):
 return ['left_minimum_second_inflection'] if any(iv==1 and qsign==-1 and psign==1 and all(sig[j]==sg for j,sg in MASK.items())for sig,iv,qsign,psign in patterns(word))else []

def check_certificates():
 data=json.loads((ROOT/'data/left_minimum_inflection_certificate.json').read_text())
 assert data['degrees']==[10,5,5]
 assert data['mask']==[[i,sg]for i,sg in sorted(MASK.items())]
 raw=raw_inequalities();bs=[];names=[]
 for row in data['polynomials']:
  name=row['name'];assert name not in names;names.append(name)
  reduced=decode(row['reduced']);factors={(0,0,0):1}
  content=row['positive_content'];assert type(content)is int and content>0
  for f in row['factors']:
   assert f['kind']=='positive_on_open_box' and type(f['exponent'])is int and f['exponent']>0
   p=decode(f['terms']);_,bb=bernstein(p)
   assert min(bb)>=0 and max(bb)>0
   factors=mul(factors,power(p,f['exponent']))
  assert scale(mul(reduced,factors),content)==raw[name]
  shape,b=bernstein(reduced,data['degrees']);bs.append(b)
 assert set(names)=={'Q0','Q1','Q3','Q6','R3','H03','H63','H67'}
 stats={'nodes':0,'leaves':0,'maxdepth':0}
 def walk(tree,arrays,depth):
  stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
  if type(tree)is int:
   assert 0<=tree<len(arrays)
   assert max(arrays[tree])<=0,('nonnegative constraint not contradicted',names[tree],depth,max(arrays[tree]))
   stats['leaves']+=1;return
  assert type(tree)is list and len(tree)==3 and tree[0]in(0,1,2)
  pairs=[split(a,shape,tree[0])for a in arrays]
  walk(tree[1],[a for a,b in pairs],depth+1)
  walk(tree[2],[b for a,b in pairs],depth+1)
 walk(data['tree'],bs,0)
 assert stats==data['statistics']
 print('EXACTLY VERIFIED coupled left-minimum/inflection obstruction:',len(names),
       'identities,',stats['leaves'],'three-dimensional Bernstein boxes.',flush=True)
 return stats

if __name__=='__main__':
 check_certificates()
 known=json.loads((ROOT/'data/certificates_n6.json').read_text())
 assert not any(reasons(w)for w in known)
 print('EXACTLY VERIFIED: no known hyperbolic witness is excluded.',flush=True)
