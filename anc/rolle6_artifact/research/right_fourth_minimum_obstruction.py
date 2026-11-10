"""Exact right-P4 / first-minimum obstruction; standard library only."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from math import comb
import json
from bernstein3 import bernstein,split,leaf
from inflection_obstruction import patterns
from triple_inflection_obstruction import identity_check as geometry_check
ROOT=Path(__file__).resolve().parents[1]
ONE={(0,0,0):F(1)};U={(1,0,0):F(1)};K={(0,1,0):F(1)};Y={(0,0,1):F(1)}

def add(*args):
 out={}
 for p in args:
  for i,c in p.items():out[i]=out.get(i,F(0))+c
 return {i:c for i,c in out.items() if c}
def scale(p,c):return {i:a*c for i,a in p.items()if a*c}
def mul(p,q):
 out={}
 for i,a in p.items():
  for j,b in q.items():
   e=tuple(x+y for x,y in zip(i,j));out[e]=out.get(e,F(0))+a*b
 return {i:c for i,c in out.items()if c}
def power(p,n):
 out=ONE
 for _ in range(n):out=mul(out,p)
 return out

def translate(p,shifts):
 out={}
 for ex,c in p.items():
  for ij in product(*(range(e+1)for e in ex)):
   v=c
   for e,i,h in zip(ex,ij,shifts):v*=comb(e,i)*h**(e-i)
   out[ij]=out.get(ij,F(0))+v
 return {i:c for i,c in out.items()if c}

def direct_minimum(p,box):
 # Independent closed formula: no use of tensor conversion or subdivision.
 degrees=tuple(max(e[i]for e in p)for i in range(3))
 basis={}
 for axis,n in enumerate(degrees):
  l,h=box[axis]
  for j in range(n+1):
   for e in range(n+1):
    basis[axis,j,e]=sum(F(comb(e,i)*comb(j,i),comb(n,i))*l**(e-i)*(h-l)**i
                         for i in range(min(e,j)+1))
 out=[]
 for j in product(*(range(n+1)for n in degrees)):
  val=0
  for ex,c in p.items():
   v=c
   for axis in range(3):v*=basis[axis,j[axis],ex[axis]]
   val+=v
  out.append(val)
 return min(out)

def check_certificates(path=None):
 geometry_check()
 data=json.loads(Path(path or ROOT/'data/right_fourth_minimum_certificate.json').read_text())
 ps=data['polynomials'];polys={name:{tuple(ij):F(c)for ij,c in ps[name]}for name in ['S','N','Q']}
 S,N,Q=map(polys.get,['S','N','Q'])
 D=add(U,scale(K,3),scale(ONE,2));V=add(ONE,K)
 assert N==add(scale(power(K,2),6),scale(K,6),scale(power(U,2),-1),scale(U,-1),ONE)
 A=add(scale(D,2),N);TN=add(mul(add(ONE,U),D),N)
 HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
 assert HN==scale(mul(V,add(scale(mul(D,V),2),scale(A,-1))),3)
 coeff={6:scale(D,2),5:scale(A,-3),4:scale(HN,-5),3:scale(mul(U,TN),10)}
 def E(x):return add(*(mul(co,power(x,n))for n,co in coeff.items()))
 EV=add(*(mul(co,power(V,n-1))for n,co in coeff.items()))
 assert mul(Q,mul(Y,add(Y,V)))==add(E(scale(Y,-1)),mul(Y,EV))
 positive=add(scale(power(K,2),3),scale(mul(K,U),2),scale(K,4),U,ONE)
 left=add(mul(power(D,5),EV),scale(add(*(mul(mul(co,power(TN,n-1)),power(D,6-n))for n,co in coeff.items())),-1))
 assert left==mul(mul(D,positive),S)
 # H>0 and r>0 imply u>3/5. r>0, u>=4/5 imply k>1/15.
 assert F(3,5)**2+F(3,5)-1==F(-1,25)
 assert 6*F(1,15)**2+6*F(1,15)-F(4,5)**2-F(4,5)+1==F(-1,75)
 for shifts in [(F(3,5),F(1,14),F(0)),(F(4,5),F(1,15),F(0))]:
  tr=translate(S,shifts)
  assert all(c>0 for c in tr.values()) and tr.get((0,0,0),0)>0
 box=[(F(3,5),F(4,5)),(F(0),F(1,14)),(F(3,5),F(2))]
 cover=data['cover'];assert cover['box']==[[str(a),str(b)]for a,b in box]
 assert cover['names']==['S','-N','Q']
 source=[S,scale(N,-1),Q];bs=[bernstein(p,box)for p in source]
 stats={'nodes':0,'leaves':0,'max_depth':0,'types':[0,0,0]};direct_bounds=[]
 def walk(tree,bs,box,depth):
  stats['nodes']+=1;stats['max_depth']=max(stats['max_depth'],depth)
  if type(tree)is int:
   assert 0<=tree<3
   m=min(bs[tree][1].values());assert (m>=0 if tree<2 else m>0)
   independent=direct_minimum(source[tree],box);assert independent==m
   stats['leaves']+=1;stats['types'][tree]+=1
   direct_bounds.append((tree,str(m)));return
  assert type(tree)is list and len(tree)==3 and type(tree[0])is int and tree[0]in range(3)
  axis=tree[0];pairs=[split(b,axis)for b in bs];mid=sum(box[axis])/2
  for side in range(2):
   b=box[:];b[axis]=(box[axis][0],mid)if side==0 else(mid,box[axis][1])
   walk(tree[1+side],[p[side]for p in pairs],b,depth+1)
 walk(cover['tree'],bs,box,0);assert stats==cover['statistics']
 # For y>=2, every coefficient in powers of (y-2) is positive on the u,k box.
 tail=translate(Q,(F(0),F(0),F(2)));tail_bounds=[]
 for j in range(5):
  p={(i,k,0):v for (i,k,l),v in tail.items()if l==j}
  b=[box[0],box[1],(F(0),F(1))];m=direct_minimum(p,b);assert m>0;tail_bounds.append(str(m))
 print('EXACTLY VERIFIED right-P4/first-minimum identities; 17 independently checked boxes; 5 positive tail coefficients.',flush=True)
 return {'cover':stats,'leaf_bounds':direct_bounds,'tail_bounds':tail_bounds}

def reasons(w):
 ps=patterns(w);sig=ps[0][0]
 if sig[5]==-1 and sig[12]==-1 and ps[1][3]==1 and ps[3][3]==1:
  return ['right_fourth_first_minimum']
 return []

if __name__=='__main__':
 check_certificates()
 known=json.loads((ROOT/'data/certificates_n6.json').read_text());unknown=json.loads((ROOT/'data/unresolved_n6.json').read_text())
 assert not any(reasons(w)for w in known)
 print('newly matched:',[w for w in unknown if reasons(w)])
