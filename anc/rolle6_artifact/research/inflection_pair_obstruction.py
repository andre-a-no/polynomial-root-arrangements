"""Exact necessary conditions at two zeros of P'', with Hermite consequences.

Polynomials use integer coefficients in k,s,t. Equality E=0 is retained.
"""
import json,math,itertools
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from left_minimum_inflection_obstruction import add,scale,mul,power,decode,bernstein,split
from inflection_obstruction import patterns as single_patterns
ROOT=Path(__file__).resolve().parents[1]

def patterns(word):
 ps=single_patterns(word)
 return [(ps[i][0],*ps[i][1:],*ps[j][1:])for i,j in combinations(range(4),2)]

def key(p):return (tuple(p[0]),*p[1:])

def raw_system(pat):
 sig,iv,uq,up,jv,vq,vp=key(pat)
 assert len(sig)==18 and all(v in(-1,1)for v in sig)
 assert 0<=iv<jv<=6 and all(v in(-1,1)for v in(uq,up,vq,vp))
 one={(0,0,0):1};K={(1,0,0):1};S={(0,1,0):1};T={(0,0,1):1}
 D=add(power(K,2),scale(K,-1),one)
 AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
 CN=add(scale(power(K,2),2),scale(K,-2),scale(one,-1));DN=scale(add(AN,CN),-1)
 BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
 ys=[AN,scale(D,-1),CN,{},D,DN];ends=[scale(D,-3)]+ys+[scale(D,3)]
 U=add(ends[iv],mul(S,add(ends[iv+1],scale(ends[iv],-1))))
 V=add(ends[jv],mul(T,add(ends[jv+1],scale(ends[jv],-1))))
 E=add(scale(BN,4),scale(mul(add(U,V),add(scale(power(D,2),6),scale(power(U,2),-1),scale(power(V,2),-1))),-1))
 A4=mul(mul(U,V),add(power(U,2),mul(U,V),power(V,2),scale(power(D,2),-6)))
 def Q(Y):return add(power(Y,5),scale(mul(power(D,2),power(Y,3)),-10),scale(mul(BN,power(Y,2)),10),scale(mul(A4,Y),5))
 def R(Y):return add(power(Y,4),scale(mul(power(D,2),power(Y,2)),-6),scale(mul(BN,Y),4),A4)
 # Off the constraint surface the residual is exactly Y*E.
 assert R(U)==mul(U,E) and R(V)==mul(V,E)
 ys.extend([U,V]);ps=list(sig[::3])+[up,vp];qs=list(sig[1::3])+[uq,vq]
 bounds={};lower=[];upper=[];meta={}
 for i,Y in enumerate(ys):
  name='W'+str(i);bounds[name]=scale(Q(Y),-6)
  (lower if qs[i]>0 else upper).append(name);meta[name]=[(3*i+1,qs[i])]if i<6 else []
 rank=[i+(i>=iv)+(i>=jv)for i in range(6)]+[iv,jv+1]
 for i in range(8):
  for j in range(i+1,8):
   if ps[i]==ps[j]:continue
   Y,Z=ys[i],ys[j]
   sums=lambda n:add(*(mul(power(Y,n-h),power(Z,h))for h in range(n+1)))
   div=add(sums(5),scale(mul(power(D,2),sums(3)),-15),scale(mul(BN,sums(2)),20),scale(mul(A4,add(Y,Z)),15))
   name='S'+str(i)+str(j);bounds[name]=scale(div,-1)
   (lower if ps[i]*(rank[i]-rank[j])>0 else upper).append(name)
   meta[name]=[(3*z,ps[z])for z in (i,j)if z<6]
 out={};masks={}
 for lo in lower:
  for hi in upper:
   name=hi+'-'+lo;out[name]=add(bounds[hi],scale(bounds[lo],-1));masks[name]=meta[hi]+meta[lo]
 for i,sg in enumerate(sig[2::3]):out['R'+str(i)]=scale(R(ys[i]),sg);masks['R'+str(i)]=[(3*i+2,sg)]
 # G=5-4A-2b^2; denominator D^6 is positive.
 out['hermite3']=add(scale(power(D,6),5),scale(mul(A4,power(D,2)),-1),scale(power(BN,2),-2));masks['hermite3']=[]
 # F(x)>0 is necessary when P(x)>0, from Delta4 and elimination of B,C.
 for i,Y in enumerate(ys):
  if ps[i]<0:continue
  f=add(scale(power(A4,2),-100),scale(mul(mul(A4,BN),Y),-120),
   scale(mul(mul(A4,power(Y,2)),power(D,2)),57),scale(mul(A4,power(D,4)),375),
   scale(mul(power(BN,2),power(Y,2)),-36),scale(mul(power(BN,2),power(D,2)),500),
   scale(mul(mul(BN,power(Y,3)),power(D,2)),100),scale(mul(mul(BN,Y),power(D,4)),-300),
   scale(mul(power(Y,6),power(D,2)),5),scale(mul(power(Y,4),power(D,4)),-75),
   scale(mul(power(Y,2),power(D,6)),90))
  out['hermite4_'+str(i)]=f;masks['hermite4_'+str(i)]=[(3*i,1)]if i<6 else []
 out['physical']=scale(CN,-1);masks['physical']=[]
 return out,E,masks


def required_mask(cert):
 sig=key(cert['pattern'])[0];mask={}
 for row in cert['polynomials']:
  name=row['name'];indices=[]
  if name in ('E+','E-','physical','hermite3'):continue
  if name.startswith('hermite4_'):
   i=int(name.split('_')[1]);indices=[3*i]if i<6 else []
  elif name.startswith('R'):indices=[3*int(name[1:])+2]
  else:
   for part in name.split('-'):
    if part[0]=='W':
     i=int(part[1:]);indices.extend([3*i+1]if i<6 else [])
    else:
     assert part[0]=='S';indices.extend(3*int(i)for i in part[1:]if int(i)<6)
  for i in indices:mask[i]=sig[i]
 return tuple(sorted(mask.items()))

def check_one(cert):
 pat=key(cert['pattern']);raw,e,masks=raw_system(pat);raw.update({'E+':e,'E-':scale(e,-1)});degrees=cert['degrees'];bs=[];shapes=[];names=[]
 assert cert['mask']==[list(x)for x in required_mask(cert)]
 assert len(degrees)==3 and all(type(n)is int and 0<=n<=20 for n in degrees)
 for row in cert['polynomials']:
  name=row['name'];assert name not in names;names.append(name);reduced=decode(row['reduced'])
  product={(0,0,0):1};content=row['positive_content'];assert type(content)is int and content>0
  for f in row['factors']:
   fp=decode(f['terms']);n=f['exponent'];assert type(n)is int and n>0
   if f['kind']=='physical':assert fp==raw['physical']
   else:
    assert f['kind']=='positive_on_open_box';_,b=bernstein(fp)
    assert min(b)>=0 and max(b)>0
   product=mul(product,power(fp,n))
  assert scale(mul(reduced,product),content)==raw[name]
  shape,b=bernstein(reduced);assert all(n-1<=d for n,d in zip(shape,degrees));bs.append(b);shapes.append(shape)
 stats={'nodes':0,'leaves':0,'maxdepth':0}
 needs={}
 def needed(tree):
  if type(tree)is int:return {tree}
  if type(tree)is dict:return {tree['i'],tree['e']}
  if id(tree)in needs:return needs[id(tree)]
  v=needed(tree[1])|needed(tree[2]);needs[id(tree)]=v;return v
 needed(cert['tree'])
 def elevate(a,shape):
  for axis,N in enumerate(degrees):
   n=shape[axis]-1
   if n==N:continue
   oldstride=(shape[1]*shape[2],shape[2],1);ns=list(shape);ns[axis]=N+1;ns=tuple(ns)
   newstride=(ns[1]*ns[2],ns[2],1);out=[F(0)]*math.prod(ns)
   other=[j for j in range(3)if j!=axis]
   weights=[[(j,F(math.comb(n,j)*math.comb(N-n,i-j),math.comb(N,i)))
             for j in range(max(0,i-(N-n)),min(n,i)+1)]for i in range(N+1)]
   for ix in itertools.product(*(range(shape[j])for j in other)):
    oldbase=sum(i*oldstride[j]for i,j in zip(ix,other));newbase=sum(i*newstride[j]for i,j in zip(ix,other))
    row=[a[oldbase+j*oldstride[axis]]for j in range(n+1)]
    for i,ws in enumerate(weights):out[newbase+i*newstride[axis]]=sum((row[j]*v for j,v in ws),F(0))
   a=out;shape=ns
  return a
 def walk(tree,arrays,depth):
  stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
  if type(tree)is dict:
   assert set(tree)=={'i','e','lambda'} and names[tree['e']]=='E+' and names[tree['i']]not in('E+','E-')
   q=F(*tree['lambda']);a=elevate(arrays[tree['i']],shapes[tree['i']]);b=elevate(arrays[tree['e']],shapes[tree['e']])
   assert max(x+q*y for x,y in zip(a,b))<=0
   stats['leaves']+=1;return
  if type(tree)is int:
   assert 0<=tree<len(bs)and tree in arrays
   a=arrays[tree]
   strict=names[tree]in('E+','E-')
   if (max(a)>=0 if strict else max(a)>0):a=elevate(a,shapes[tree])
   assert (max(a)<0 if strict else max(a)<=0),(pat,names[tree],depth,max(a))
   stats['leaves']+=1;return
  assert type(tree)is list and len(tree)==3 and tree[0]in(0,1,2)
  ln=needed(tree[1])
  rn=needed(tree[2])
  pairs={i:split(arrays[i],shapes[i],tree[0])for i in ln|rn}
  walk(tree[1],{i:pairs[i][0]for i in ln},depth+1)
  walk(tree[2],{i:pairs[i][1]for i in rn},depth+1)
 walk(cert['tree'],dict(enumerate(bs)),0);assert stats==cert['statistics']
 return stats

def check_certificates():
 from hermite_criterion import check_identities
 check_identities()
 data=json.loads((ROOT/'data/inflection_pair_certificates.json').read_text());seen=set();leaves=0
 for cert in data['certificates']:
  p=key(cert['pattern']);assert p not in seen;seen.add(p);leaves+=check_one(cert)['leaves']
 print('EXACTLY VERIFIED',len(seen),'two-inflection certificates and',leaves,'three-dimensional Bernstein boxes on the equality surface.',flush=True)
 return seen

_certified=None
def reasons(word):
 global _certified
 if _certified is None:
  cs=json.loads((ROOT/'data/inflection_pair_certificates.json').read_text())['certificates']
  _certified=[(key(c['pattern'])[1:],required_mask(c),i)for i,c in enumerate(cs)]
 return ['inflection_pair_'+str(i)for p in patterns(word)for tail,mask,i in _certified
         if p[1:]==tail and all(p[0][j]==sg for j,sg in mask)]

if __name__=='__main__':
 check_certificates()
 known=json.loads((ROOT/'data/certificates_n6.json').read_text())
 assert not any(reasons(w)for w in known)
 print('EXACTLY VERIFIED: no known real-rooted witness excluded.',flush=True)
