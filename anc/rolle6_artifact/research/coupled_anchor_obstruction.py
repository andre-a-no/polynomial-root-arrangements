"""Coupled zero of P'' and zero of P': three-dimensional exact obstructions."""
import json,math,itertools
from fractions import Fraction as F
from pathlib import Path
from inflection_obstruction import patterns as inflection_patterns
from stationary_obstruction import patterns as critical_patterns
from left_minimum_inflection_obstruction import add,scale,mul,power,decode,bernstein,split
ROOT=Path(__file__).resolve().parents[1]

def patterns(word):
 pos=[[i for i,c in enumerate(word)if int(c)==j]for j in range(6)]
 if pos[3][1]>pos[5][0]:return patterns(word[::-1])
 return [(ip[0],ip[1],ip[2],ip[3],cp[1],cp[2],int(pos[2][i]<pos[1][j]))
         for i,ip in enumerate(inflection_patterns(word))for j,cp in enumerate(critical_patterns(word))]

def key(p):return(tuple(p[0]),*p[1:])

def required_mask(cert):
 pat=key(cert['pattern']);sig=pat[0];mask={}
 for row in cert['polynomials']:
  name=row['name']
  if name=='physical':continue
  if name[0]in'QR':
   j=int(name[1:])
   if j<6:mask[3*j+(1 if name[0]=='Q'else 2)]=sig[3*j+(1 if name[0]=='Q'else 2)]
  else:
   assert name[0]=='H';i,j=map(int,name[1:])
   for index,sg in [(i,1),(j,-1)]:
    if index<6:
     assert sig[3*index]==sg;mask[3*index]=sg
    elif index==6:assert pat[3]==sg
    else:assert index==7 and pat[5]==sg
 return tuple(sorted(mask.items()))

def raw_inequalities(pat):
 sig,wi,wq,wp,ui,up,order=key(pat)
 assert len(sig)==18 and all(v in(-1,1)for v in sig)
 assert wi in range(7)and ui in range(7)and all(v in(-1,1)for v in(wq,wp,up))and order in(0,1)
 if wi!=ui:assert order==int(wi<ui)
 one={(0,0,0):1};K={(1,0,0):1};S={(0,1,0):1};T={(0,0,1):1}
 D=add(power(K,2),scale(K,-1),one)
 AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
 CN=add(scale(power(K,2),2),scale(K,-2),scale(one,-1));DN=scale(add(AN,CN),-1)
 BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
 assert add(mul(AN,CN),mul(AN,DN),mul(CN,DN),scale(power(D,2),3))=={}
 assert add(mul(mul(AN,CN),DN),BN)=={}
 ys=[AN,scale(D,-1),CN,{},D,DN]
 wb=[scale(D,-3)]+ys+[scale(D,3)];ub=[scale(D,-4)]+ys+[scale(D,4)]
 W=add(wb[wi],mul(S,add(wb[wi+1],scale(wb[wi],-1))))
 lo,hi=ub[ui:ui+2]
 if wi==ui:
  if order:lo=W
  else:hi=W
 U=add(lo,mul(T,add(hi,scale(lo,-1))))
 A4=add(scale(power(W,4),-1),scale(mul(power(D,2),power(W,2)),6),scale(mul(BN,W),-4))
 B5=add(scale(power(U,5),-1),scale(mul(power(D,2),power(U,3)),10),scale(mul(BN,power(U,2)),-10),scale(mul(A4,U),-5))
 def P(Y):return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20),scale(mul(A4,power(Y,2)),15),scale(mul(B5,Y),6))
 def Q(Y):return add(power(Y,5),scale(mul(power(D,2),power(Y,3)),-10),scale(mul(BN,power(Y,2)),10),scale(mul(A4,Y),5),B5)
 def R(Y):return add(power(Y,4),scale(mul(power(D,2),power(Y,2)),-6),scale(mul(BN,Y),4),A4)
 assert R(W)=={} and Q(U)=={}
 ys.extend([W,U]);out={}
 for j,sg in enumerate(list(sig[1::3])+[wq]):out[f'Q{j}']=scale(Q(ys[j]),sg)
 for j,sg in enumerate(sig[2::3]):out[f'R{j}']=scale(R(ys[j]),sg)
 out['R7']=scale(R(U),-up)
 ps=list(sig[::3])+[wp,up]
 for i in range(8):
  for j in range(8):
   if ps[i]>0 and ps[j]<0:out[f'H{i}{j}']=add(P(ys[i]),scale(P(ys[j]),-1))
 out['physical']=scale(CN,-1)
 return out

def check_one(cert):
 pat=key(cert['pattern']);raw=raw_inequalities(pat);degrees=cert['degrees'];bs=[];shapes=[];names=[]
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
  if type(tree)is int:
   assert 0<=tree<len(bs)and tree in arrays
   a=arrays[tree]
   if max(a)>0:a=elevate(a,shapes[tree])
   assert max(a)<=0,(pat,names[tree],depth,max(a))
   stats['leaves']+=1;return
  assert type(tree)is list and len(tree)==3 and tree[0]in(0,1,2)
  ln={tree[1]}if type(tree[1])is int else needs[id(tree[1])]
  rn={tree[2]}if type(tree[2])is int else needs[id(tree[2])]
  pairs={i:split(arrays[i],shapes[i],tree[0])for i in ln|rn}
  walk(tree[1],{i:pairs[i][0]for i in ln},depth+1)
  walk(tree[2],{i:pairs[i][1]for i in rn},depth+1)
 walk(cert['tree'],dict(enumerate(bs)),0);assert stats==cert['statistics']
 return stats

def check_certificates():
 data=json.loads((ROOT/'data/coupled_anchor_certificates.json').read_text());seen=set();leaves=0
 for c in data['certificates']:
  p=key(c['pattern']);assert p not in seen;seen.add(p);leaves+=check_one(c)['leaves']
 print('EXACTLY VERIFIED',len(seen),'coupled inflection/stationary patterns and',leaves,'three-dimensional Bernstein boxes.',flush=True)
 return seen

_certified=None
def reasons(word):
 global _certified
 if _certified is None:
  data=json.loads((ROOT/'data/coupled_anchor_certificates.json').read_text())
  _certified=[(key(c['pattern'])[1:],required_mask(c),i)for i,c in enumerate(data['certificates'])]
 return ['coupled_anchor_'+str(i)for p in patterns(word)for tail,mask,i in _certified
         if p[1:]==tail and all(p[0][j]==sg for j,sg in mask)]

if __name__=='__main__':
 check_certificates()
 known=json.loads((ROOT/'data/certificates_n6.json').read_text())
 assert not any(reasons(w)for w in known)
 print('EXACTLY VERIFIED: no known hyperbolic witness excluded.',flush=True)
