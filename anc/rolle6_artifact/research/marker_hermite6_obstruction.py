"""Exact full-H6 directions with a necessary B-chord bound in two parameters.
The chord relaxation is only a necessary condition. Certificates prove
infeasibility; failed searches and positive relaxed points prove no word status.
"""
import json
from fractions import Fraction as F
from math import comb,lcm
from pathlib import Path
from two_parameter_obstruction import add,scale,mul,power,bernstein2,split2
from stationary_obstruction import dense
from inflection_obstruction import patterns
from marker_hermite_obstruction import tables,BOX,decode,checked_reduction
ROOT=Path(__file__).resolve().parents[1]

def b_tables(sig,rows,D,BN,A4):
 bounds={f'Q{j}':(scale(rows[f'Q{j}'][2],-sig[3*j+1]),sig[3*j+1])for j in range(6)}
 bounds.update(Blo=(scale(power(D,5),-192),1),Bhi=(scale(power(D,5),192),-1))
 ys=[scale(rows[f'P{j}'][0],sig[3*j])for j in range(6)]
 for i,Y in enumerate(ys):
  for j,Z in enumerate(ys[i+1:],i+1):
   if sig[3*i]==sig[3*j]:continue
   sums=lambda n:add(*(mul(power(Y,n-h),power(Z,h))for h in range(n+1)))
   bd=scale(add(sums(5),scale(mul(power(D,2),sums(3)),-15),scale(mul(BN,sums(2)),20),scale(mul(A4,add(Y,Z)),15)),-1)
   direction=-sig[3*i];name=f'P{i}{j}';bounds[name]=(bd,direction)
   rows[name]=({(0,0):direction},{},scale(bd,-direction))
 return bounds

def univariate_division(p,divisor,axis):
 # The divisor depends on one variable only; retain exact rational coefficients.
 q={};r={};other=1-axis
 ds={m[axis]:F(v)for m,v in divisor.items()};assert all(m[other]==0 for m in divisor)
 for h in {m[other]for m in p}:
  rem={m[axis]:F(v)for m,v in p.items()if m[other]==h};quo={}
  while rem and max(rem)>=max(ds):
   e=max(rem)-max(ds);v=rem[max(rem)]/ds[max(ds)];quo[e]=quo.get(e,F(0))+v
   for i,z in ds.items():
    rem[i+e]=rem.get(i+e,F(0))-v*z
    if not rem[i+e]:del rem[i+e]
  for target,source in((q,quo),(r,rem)):
   for i,v in source.items():
    if v:target[(i,h)if axis==0 else(h,i)]=v
 assert add(mul(q,divisor),r)==p
 return q,r

def pair_names(pair_index,tests):
 i,j=tests
 return [f'test_{h}_l{s}'for h in(i,j)for s in(1,2)]+[f'pair_{pair_index}_{s}'for s in('w1','w2','value')]

def raw_system(sig,tests,chart='quartic',pairs=()):
 raw,rows,D,BN,_=tables(sig,chart);A4=scale(raw['R3'],sig[11]);bounds=b_tables(sig,rows,D,BN,A4)
 one={(0,0):1};co=[one,{},scale(power(D,2),-15),scale(BN,20),scale(A4,15),{},{}];mom=[scale(one,6)]
 for n in range(1,11):
  z=[mul(co[j],mom[n-j])for j in range(1,min(n-1,6)+1)]
  if n<=6:z.append(scale(co[n],n))
  mom.append(scale(add(*z),-1))
 mom=[mul(v,power(D,10-n))for n,v in enumerate(mom)];determinants={}
 for j,test in enumerate(tests):
  v=test['vector'];assert len(v)==6 and all(type(z)is int for z in v)and any(v)
  x,y=test['rows'];assert x in rows and y in rows and x!=y
  sg=test['det_sign'];assert sg in(-1,1)
  lo,sl=bounds[test['lower_bound']];hi,sh=bounds[test['upper_bound']];assert sl==1 and sh==-1
  v0,v1,v2,v3,v4,v5=v;K=72000*v5*v5
  a0=-1200*(v0*v5+v1*v4+v2*v3)-25200*(v2*v5+v3*v4)-486000*v4*v5
  ab=19200*(2*v3*v5+v4*v4)+720000*v5*v5
  c0=-720*(2*v1*v5+2*v2*v4+v3*v3)-14400*(2*v3*v5+v4*v4)-270000*v5*v5
  a4=add(scale(power(D,4),a0),scale(mul(BN,D),ab),scale(A4,32400*v4*v5))
  c4=add(scale(power(D,4),c0),scale(mul(BN,D),43200*v4*v5),scale(A4,18000*v5*v5))
  a=add(scale(mul(D,a4),120),scale(add(lo,hi),K));c=scale(c4,120)
  d10=add(*(scale(mom[h+i],v[h]*v[i])for h in range(6)for i in range(6)))
  d=add(scale(d10,14400),scale(mul(lo,hi),-K))
  # Bound the affine form at BOTH endpoints of the B interval. This
  # avoids a case split on a-cY, whose zero curve is hard to tile.
  for name in [f'P{i}'for i in range(6)]+['Clo','Chi']:
   rb,rc,rd=rows[name];assert rc in({(0,0):1},{(0,0):-1})
   direction=rc[(0,0)];Y=scale(rb,direction);Z=scale(rd,direction)
   aa=add(a,scale(mul(c,Y),-1));dd=add(d,scale(mul(c,Z),-1))
   raw[f'test_{j}_ep_{name}_sign']=scale(c,direction)
   raw[f'test_{j}_ep_{name}_lo']=add(dd,mul(aa,lo))
   raw[f'test_{j}_ep_{name}_hi']=add(dd,mul(aa,hi))
  b1,c1,d1=rows[x];b2,c2,d2=rows[y]
  det=scale(add(mul(b1,c2),scale(mul(b2,c1),-1)),sg)
  lam1=scale(add(scale(mul(a,c2),-1),mul(c,b2)),sg)
  lam2=scale(add(scale(mul(c,b1),-1),mul(a,c1)),sg)
  assert add(mul(a,det),mul(lam1,b1),mul(lam2,b2))=={}
  assert add(mul(c,det),mul(lam1,c1),mul(lam2,c2))=={}
  raw[f'test_{j}_l1']=scale(lam1,-1);raw[f'test_{j}_l2']=scale(lam2,-1)
  raw[f'test_{j}_value']=add(mul(d,det),mul(lam1,d1),mul(lam2,d2))
  n=raw[f'test_{j}_value'];p0={i:F(v)for (i,h),v in n.items()if h==0}
  div={i:F(v)for (i,h),v in raw['physical'].items()};rem=dict(p0);quo={}
  while rem and max(rem)>=max(div):
   e=max(rem)-max(div);v=rem[max(rem)]/div[max(div)];quo[e]=quo.get(e,F(0))+v
   for h,w in div.items():
    rem[h+e]=rem.get(h+e,F(0))-v*w
    if not rem[h+e]:del rem[h+e]
  if not rem:
   den=lcm(*(v.denominator for v in quo.values()))if quo else 1
   q={(i,0):int(v*den)for i,v in quo.items()if v};r={(i,h-1):v for (i,h),v in n.items()if h>0}
   assert add(scale(mul(raw['physical'],q),F(1,den)),mul({(0,1):1},r))==n
   raw[f'test_{j}_value0']=q;raw[f'test_{j}_value1']=r
  determinants[f'det_{j}']=det
 for h,pair in enumerate(pairs):
  assert len(pair)==2 and all(type(j)is int and 0<=j<len(tests)for j in pair)and pair[0]!=pair[1]
  qs=[];fs=[];ns=[];T={(0,1):15,(0,0):-7}
  for j in pair:
   n=raw[f'test_{j}_value'];q,r=univariate_division(n,raw['physical'],0);f,rem=univariate_division(r,T,1)
   assert not rem
   den=lcm(*(v.denominator for v in list(q.values())+list(f.values())))
   q={m:int(v*den)for m,v in q.items()};f={m:int(v*den)for m,v in f.items()};n=scale(n,den)
   assert add(mul(raw['physical'],q),mul(T,f))==n
   qs.append(q);fs.append(f);ns.append(n)
  w1=scale(fs[1],-1);w2=fs[0];value=add(mul(w1,qs[0]),mul(w2,qs[1]))
  assert add(mul(w1,ns[0]),mul(w2,ns[1]))==mul(raw['physical'],value)
  raw[f'pair_{h}_w1']=scale(w1,-1);raw[f'pair_{h}_w2']=scale(w2,-1);raw[f'pair_{h}_value']=value
 return raw,determinants

def clause_names(sig,j,raw):
 out=[[f'test_{j}_{s}'for s in('l1','l2','value')]]
 if f'test_{j}_value0'in raw:out.append([f'test_{j}_{s}'for s in('l1','l2','value0','value1')])
 cb=[(f'P{i}',sig[3*i])for i in range(6)]+[('Clo',1),('Chi',-1)]
 for name,sg in cb:out.append([f'test_{j}_ep_{name}_{s}'for s in('sign','lo','hi')])
 for x,sx in cb:
  for y,sy in cb:
   if sx==1 and sy==-1:
    out.append([f'test_{j}_ep_{name}_{s}'for name in(x,y)for s in('lo','hi')])
 return out

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
  if type(t)is dict:
   if 'test'in t:used_tests.add(t['test'])
   else:used_tests.update(cert['pairs'][t['pair']])
  elif type(t)is list:walk(t[1]);walk(t[2])
 walk(cert['tree'])
 for row in cert['polynomials']:
  name=row['name']
  if name.startswith('R'):add_index(3*int(name[1:])+2)
  elif name.startswith('bound_'):
   for part in name.split('_')[1:]:bound(part)
  elif '_ep_'in name:bound(name.split('_ep_')[1].rsplit('_',1)[0])
 for j in used_tests:
  for name in cert['tests'][j]['rows']:bound(name)
  bound(cert['tests'][j]['lower_bound']);bound(cert['tests'][j]['upper_bound'])
 return [[i,v]for i,v in sorted(mask.items())]

def check_one(cert):
 assert cert['mask']==required_mask(cert)
 sig=tuple(cert['signature']);raw,dets=raw_system(sig,cert['tests'],cert.get('chart','quartic'),cert.get('pairs',[]))
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
   assert 0<=t<len(names)and not names[t].startswith(('test_','pair_'))
   assert max(v for row in arr[t]for v in row)<=0
   stats['leaves']+=1;return
  if type(t)is dict:
   if 'pair'in t:
    assert set(t)=={'pair','indices'};j=t['pair'];assert type(j)is int and 0<=j<len(cert['pairs'])
    assert [names[i]for i in t['indices']]==pair_names(j,cert['pairs'][j])
    assert all(max(v for row in arr[i]for v in row)<0 for i in t['indices'][4:6])
   else:
    assert set(t)=={'test','indices'};j=t['test'];assert type(j)is int and 0<=j<len(cert['tests'])
    assert [names[i]for i in t['indices']]in clause_names(sig,j,raw)
   assert all(max(v for row in arr[i]for v in row)<=0 for i in t['indices'])
   stats['leaves']+=1;return
  assert type(t)is list and len(t)==3 and t[0]in(0,1)
  ln,rn=need(t[1]),need(t[2]);pairs={i:split2(arr[i],t[0])for i in ln|rn}
  walk(t[1],{i:pairs[i][0]for i in ln},depth+1);walk(t[2],{i:pairs[i][1]for i in rn},depth+1)
 walk(cert['tree'],dict(enumerate(bs)),0);assert stats==cert['statistics']
 return stats

def check_certificates():
 from hermite_fiber_obstruction import check_full_identities as check_identities
 check_identities();data=json.loads((ROOT/'data/marker_hermite6_certificates.json').read_text());seen=set();n=0
 for cert in data['certificates']:
  sig=tuple(cert['signature']);assert sig not in seen;seen.add(sig);n+=check_one(cert)['leaves']
 print('EXACTLY VERIFIED',len(seen),'full-H6-chord marker signatures and',n,'two-dimensional Bernstein boxes.',flush=True)
 return seen

_certified=None
def reasons(word):
 global _certified
 if _certified is None:_certified=[c['mask']for c in json.loads((ROOT/'data/marker_hermite6_certificates.json').read_text())['certificates']]
 sig=patterns(word)[0][0]
 return ['marker_full_H6_chord']if any(all(sig[i]==v for i,v in mask)for mask in _certified)else[]

if __name__=='__main__':check_certificates()
