"""Exact H5-vector and two-row Farkas obstructions in the (k,A) chart.
Independent rational checker; numerical programs only propose certificates.
"""
import json
from math import comb,lcm
from fractions import Fraction as F
from pathlib import Path
from two_parameter_obstruction import add,scale,mul,power,bernstein2,split2
from stationary_obstruction import dense
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
BOX=(F(1),F(11,8),F(0),F(1))

def tables(sig,chart='quartic'):
 assert len(sig)==18 and all(v in(-1,1)for v in sig)
 one={(0,0):1};K={(1,0):1};S={(0,1):1}
 D=add(power(K,2),scale(K,-1),one)
 AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
 CN=add(scale(power(K,2),2),scale(K,-2),scale(one,-1));DN=scale(add(AN,CN),-1)
 BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
 assert chart in('quartic','qa_Bpos','qa_qd_Bpos','qd_Qc')
 if chart!='quartic':
  D=scale(D,5);AN=scale(AN,5);CN=scale(CN,5);DN=scale(DN,5);BN=scale(BN,125)
 ys=[AN,scale(D,-1),CN,{},D,DN]
 def V(Y):return add(scale(power(Y,4),-1),scale(mul(power(D,2),power(Y,2)),6),scale(mul(BN,Y),-4))
 L=V(CN);U=V(AN)
 if chart!='quartic':
  if chart!='qd_Qc':assert sig[1]==-1 and sig[10]==1
  def FQ(Y):
   y4=power(Y,4);assert all(v%5==0 for v in y4.values())
   return add({m:-v//5 for m,v in y4.items()},scale(mul(power(Y,2),power(D,2)),2),scale(mul(BN,Y),-2))
  L=FQ(AN)
  if chart=='qd_Qc':
   assert sig[16]==1 and sig[7]==-1
   sums=lambda n:add(*(mul(power(DN,n-i),power(CN,i))for i in range(n+1)))
   s4=sums(4);assert all(v%5==0 for v in s4.values())
   L=add({m:-v//5 for m,v in s4.items()},scale(mul(power(D,2),sums(2)),2),scale(mul(BN,add(DN,CN)),-2))
  if chart=='qa_qd_Bpos':
   assert sig[16]==-1;U=FQ(DN)
 A4=add(L,mul(S,add(U,scale(L,-1))))
 def P(Y):return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20),scale(mul(A4,power(Y,2)),15))
 def Q(Y):return scale(add(power(Y,5),scale(mul(power(D,2),power(Y,3)),-10),scale(mul(BN,power(Y,2)),10),scale(mul(A4,Y),5)),6)
 rows={}
 for j,Y in enumerate(ys):
  ps,qs,rs=sig[3*j:3*j+3]
  rows['P'+str(j)]=(scale(Y,ps),scale(one,ps),scale(P(Y),ps))
  rows['Q'+str(j)]=(scale(one,qs),{},scale(Q(Y),qs))
 rows.update(Blo=(one,{},scale(power(D,5),192)),Bhi=(scale(one,-1),{},scale(power(D,5),192)),
             Clo=({},one,scale(power(D,6),125)),Chi=({},scale(one,-1),scale(power(D,6),125)))
 raw={'physical':scale(CN,-1),'chart_width':add(U,scale(L,-1))}
 for j,Y in enumerate(ys):raw['R'+str(j)]=scale(add(A4,scale(V(Y),-1)),sig[3*j+2])
 raw['H3']=add(scale(power(D,6),5),scale(power(BN,2),-2),scale(mul(A4,power(D,2)),-1))
 raw['moment1']=add(scale(power(D,4),250),scale(mul(BN,D),79),scale(A4,-50))
 raw['moment2']=add(scale(power(D,4),1021),scale(mul(BN,D),698),scale(A4,-125))
 # Eliminate C from all marker signs, then compare the necessary B bounds.
 bounds=[(scale(Q(Y),-1),sig[3*j+1],'Q'+str(j))for j,Y in enumerate(ys)]
 bounds.extend([(scale(power(D,5),-192),1,'Blo'),(scale(power(D,5),192),-1,'Bhi')])
 for j,Y in enumerate(ys):
  for h,Z in enumerate(ys[j+1:],j+1):
   if sig[3*j]==sig[3*h]:continue
   sums=lambda n:add(*(mul(power(Y,n-i),power(Z,i))for i in range(n+1)))
   bd=scale(add(sums(5),scale(mul(power(D,2),sums(3)),-15),scale(mul(BN,sums(2)),20),scale(mul(A4,add(Y,Z)),15)),-1)
   assert mul(bd,add(Y,scale(Z,-1)))==add(P(Z),scale(P(Y),-1))
   bounds.append((bd,-sig[3*j],'P'+str(j)+str(h)))
 for lo,sg,ln in bounds:
  if sg!=1:continue
  for hi,sg2,hn in bounds:
   if sg2==-1:raw['bound_'+ln+'_'+hn]=add(hi,scale(lo,-1))
 co=[one,{},scale(power(D,2),-15),scale(BN,20),scale(A4,15),{},{}];mom=[scale(one,6)]
 for n in range(1,9):
  z=[mul(co[j],mom[n-j])for j in range(1,min(n-1,6)+1)]
  if n<=6:z.append(scale(co[n],n))
  mom.append(scale(add(*z),-1))
 mom=[mul(v,power(D,8-n))for n,v in enumerate(mom)]
 return raw,rows,D,BN,mom

def raw_system(sig,tests,chart='quartic'):
 raw,rows,D,BN,mom=tables(sig,chart);determinants={}
 for j,test in enumerate(tests):
  v=test['vector'];assert len(v)==5 and all(type(z)is int for z in v)and any(v)
  x,y=test['rows'];assert x in rows and y in rows and x!=y
  sg=test['det_sign'];assert sg in(-1,1)
  a0=-1200*(v[1]*v[4]+v[2]*v[3])-25200*v[3]*v[4]
  a=add(scale(power(D,3),a0),scale(BN,19200*v[4]**2))
  c=scale(power(D,2),-720*(2*v[2]*v[4]+v[3]**2+20*v[4]**2))
  d=scale(add(*(scale(mom[h+i],v[h]*v[i])for h in range(5)for i in range(5))),120)
  b1,c1,d1=rows[x];b2,c2,d2=rows[y]
  det=scale(add(mul(b1,c2),scale(mul(b2,c1),-1)),sg)
  lam1=scale(add(scale(mul(a,c2),-1),mul(c,b2)),sg)
  lam2=scale(add(scale(mul(c,b1),-1),mul(a,c1)),sg)
  assert add(mul(a,det),mul(lam1,b1),mul(lam2,b2))=={}
  assert add(mul(c,det),mul(lam1,c1),mul(lam2,c2))=={}
  raw[f'test_{j}_l1']=scale(lam1,-1);raw[f'test_{j}_l2']=scale(lam2,-1)
  raw[f'test_{j}_value']=add(mul(d,det),mul(lam1,d1),mul(lam2,d2))
  # Boundary-aware decomposition n = physical*q(k) + s*r(k,s).
  # Each factor physical and s is strictly positive in the original chart.
  n=raw[f'test_{j}_value'];p0={i:F(v)for (i,h),v in n.items()if h==0}
  div={i:F(v)for (i,h),v in raw['physical'].items()};rem=dict(p0);quo={}
  while rem and max(rem)>=max(div):
   e=max(rem)-max(div);v=rem[max(rem)]/div[max(div)];quo[e]=quo.get(e,F(0))+v
   for h,w in div.items():
    rem[h+e]=rem.get(h+e,F(0))-v*w
    if not rem[h+e]:del rem[h+e]
  if not rem:
   den=lcm(*(v.denominator for v in quo.values()))if quo else 1
   q={(i,0):int(v*den)for i,v in quo.items()if v}
   r={(i,h-1):v for (i,h),v in n.items()if h>0}
   assert scale(mul(raw['physical'],q),F(1,den))=={(i,0):v for i,v in p0.items()}
   assert add(scale(mul(raw['physical'],q),F(1,den)),mul({(0,1):1},r))==n
   raw[f'test_{j}_value0']=q;raw[f'test_{j}_value1']=r
  determinants[f'det_{j}']=det
 return raw,determinants

def decode(terms):
 out={}
 for m,c in terms:
  assert len(m)==2 and all(type(i)is int and i>=0 for i in m)and type(c)is int
  assert tuple(m)not in out;out[tuple(m)]=c
 return {m:c for m,c in out.items()if c}

def checked_reduction(record,raw,physical):
 red=decode(record['reduced']);product={(0,0):1};content=record['positive_content']
 assert type(content)is int and content>0
 for factor in record['factors']:
  p=decode(factor['terms']);n=factor['exponent'];assert type(n)is int and n>0
  if factor['kind']=='physical':assert p==physical
  else:
   assert factor['kind']=='positive_on_open_box'
   b=bernstein2(dense(p),*BOX);flat=[v for row in b for v in row]
   assert min(flat)>=0 and max(flat)>0
  product=mul(product,power(p,n))
 assert scale(mul(red,product),content)==raw
 return bernstein2(dense(red),*BOX)

def required_mask(cert):
 sig=cert['signature'];mask={};chart=cert.get('chart','quartic')
 def add_index(i):mask[i]=sig[i]
 if chart in('qa_Bpos','qa_qd_Bpos'):
  add_index(1);add_index(10)
  if chart=='qa_qd_Bpos':add_index(16)
 if chart=='qd_Qc':add_index(16);add_index(7)
 def bound(name):
  if name in('Blo','Bhi','Clo','Chi'):return
  if name.startswith('Q'):add_index(3*int(name[1:])+1)
  else:
   assert name.startswith('P')
   for j in name[1:]:add_index(3*int(j))
 used_tests=set()
 def walk(t):
  if type(t)is dict:used_tests.add(t['test'])
  elif type(t)is list:walk(t[1]);walk(t[2])
 walk(cert['tree'])
 for row in cert['polynomials']:
  name=row['name']
  if name.startswith('R'):add_index(3*int(name[1:])+2)
  elif name.startswith('bound_'):
   for part in name.split('_')[1:]:bound(part)
 for j in used_tests:
  for name in cert['tests'][j]['rows']:bound(name)
 return [[i,v]for i,v in sorted(mask.items())]

def check_one(cert):
 assert cert['mask']==required_mask(cert)
 sig=tuple(cert['signature']);raw,dets=raw_system(sig,cert['tests'],cert.get('chart','quartic'))
 for r in cert['determinants']:
  b=checked_reduction(r,dets[r['name']],raw['physical']);flat=[v for row in b for v in row]
  assert min(flat)>=0 and max(flat)>0
 assert {r['name']for r in cert['determinants']}==set(dets)
 degrees=cert['degrees'];assert len(degrees)==2 and all(type(n)is int and 0<=n<=32 for n in degrees)
 def elevate(bs):
  n=len(bs)-1;m=len(bs[0])-1;N,M=degrees
  assert n<=N and m<=M
  a=[[sum((bs[j][h]*F(comb(n,j)*comb(N-n,i-j),comb(N,i))for j in range(max(0,i-(N-n)),min(n,i)+1)),F(0))for h in range(m+1)]for i in range(N+1)]
  return [[sum((a[i][j]*F(comb(m,j)*comb(M-m,h-j),comb(M,h))for j in range(max(0,h-(M-m)),min(m,h)+1)),F(0))for h in range(M+1)]for i in range(N+1)]
 names=[];bs=[]
 for r in cert['polynomials']:
  name=r['name'];assert name in raw and name not in names
  names.append(name);bs.append(elevate(checked_reduction(r,raw[name],raw['physical'])))
 stats={'nodes':0,'leaves':0,'maxdepth':0}
 needed={}
 def need(t):
  if type(t)is int:return {t}
  if type(t)is dict:return set(t['indices'])
  if id(t)not in needed:needed[id(t)]=need(t[1])|need(t[2])
  return needed[id(t)]
 def walk(t,arr,depth):
  stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
  if type(t)is int:
   assert 0<=t<len(names)and not names[t].startswith('test_')
   assert max(v for row in arr[t]for v in row)<=0
   stats['leaves']+=1;return
  if type(t)is dict:
   assert set(t)=={'test','indices'};j=t['test'];assert type(j)is int and 0<=j<len(cert['tests'])
   suffix=('l1','l2','value')if len(t['indices'])==3 else('l1','l2','value0','value1')
   assert [names[i]for i in t['indices']]==[f'test_{j}_{s}'for s in suffix]
   assert all(max(v for row in arr[i]for v in row)<=0 for i in t['indices'])
   stats['leaves']+=1;return
  assert type(t)is list and len(t)==3 and t[0]in(0,1)
  ln,rn=need(t[1]),need(t[2]);pairs={i:split2(arr[i],t[0])for i in ln|rn}
  walk(t[1],{i:pairs[i][0]for i in ln},depth+1);walk(t[2],{i:pairs[i][1]for i in rn},depth+1)
 walk(cert['tree'],dict(enumerate(bs)),0);assert stats==cert['statistics']
 return stats

def check_certificates():
 from hermite_fiber_obstruction import check_identities
 check_identities();data=json.loads((ROOT/'data/marker_hermite_certificates.json').read_text());seen=set();n=0
 for cert in data['certificates']:
  sig=tuple(cert['signature']);assert sig not in seen;seen.add(sig);n+=check_one(cert)['leaves']
 print('EXACTLY VERIFIED',len(seen),'full-H5 marker signatures and',n,'two-dimensional Bernstein boxes.',flush=True)
 return seen

_certified=None
def reasons(word):
 global _certified
 if _certified is None:_certified=[c['mask']for c in json.loads((ROOT/'data/marker_hermite_certificates.json').read_text())['certificates']]
 sig=patterns(word)[0][0]
 return ['marker_full_H5']if any(all(sig[i]==v for i,v in mask)for mask in _certified)else[]

if __name__=='__main__':check_certificates()
