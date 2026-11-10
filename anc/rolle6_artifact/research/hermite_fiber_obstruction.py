"""Hermite H5/H6 test vectors combined with polynomial sign constraints.

H5 clauses have two members; H6 chord clauses have three members.
The vectors are integers; no semidefinite numerical output is trusted.
"""
import json,math,itertools
from fractions import Fraction as F
from pathlib import Path
from left_minimum_inflection_obstruction import add,scale,mul,power,decode,bernstein,split
from inflection_pair_obstruction import raw_system as base_system,patterns,key,required_mask as base_mask
ROOT=Path(__file__).resolve().parents[1]

def extra_system(pat,tests):
 sig,iv,uq,up,jv,vq,vp=key(pat)
 one={(0,0,0):1};K={(1,0,0):1};S={(0,1,0):1};T={(0,0,1):1}
 D=add(power(K,2),scale(K,-1),one);AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
 CN=add(scale(power(K,2),2),scale(K,-2),scale(one,-1));DN=scale(add(AN,CN),-1)
 BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
 ys=[AN,scale(D,-1),CN,{},D,DN];ends=[scale(D,-3)]+ys+[scale(D,3)]
 U=add(ends[iv],mul(S,add(ends[iv+1],scale(ends[iv],-1))))
 V=add(ends[jv],mul(T,add(ends[jv+1],scale(ends[jv],-1))))
 A4=mul(mul(U,V),add(power(U,2),mul(U,V),power(V,2),scale(power(D,2),-6)))
 ys.extend([U,V]);ps=list(sig[::3])+[up,vp];qs=list(sig[1::3])+[uq,vq]
 def Q(Y):return add(power(Y,5),scale(mul(power(D,2),power(Y,3)),-10),scale(mul(BN,power(Y,2)),10),scale(mul(A4,Y),5))
 def P(Y):return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20),scale(mul(A4,power(Y,2)),15))
 # D^8 times Newton sums at B=C=0. They are reconstructed by recurrence.
 den=[one,D,power(D,2),power(D,3),power(D,4),power(D,5),power(D,6),power(D,7),power(D,8)]
 # Scaled moments M_j=D^j s_j, for the monic polynomial 120P.
 coeff=[one,{},scale(power(D,2),-15),scale(BN,20),scale(A4,15),{},{}]
 moments=[scale(one,6)]
 for n in range(1,9):
  rhs=add(*(mul(coeff[j],moments[n-j])for j in range(1,min(n-1,6)+1)))
  if n<=6:rhs=add(rhs,scale(coeff[n],n))
  moments.append(scale(rhs,-1))
 moments=[mul(m,den[8-j])for j,m in enumerate(moments)]
 out={};masks={}
 for index,test in enumerate(tests):
  vec=test['vector'];i=test['p_node'];j=test['q_node']
  assert len(vec)==5 and all(type(z)is int for z in vec)and any(vec)
  assert type(i)is int and type(j)is int and -1<=i<8 and -1<=j<8
  v0,v1,v2,v3,v4=vec
  a0=-1200*(v1*v4+v2*v3)-25200*v3*v4;a1=19200*v4*v4
  c=-720*(2*v2*v4+v3*v3+20*v4*v4)
  assert (i==-1 and c==0)or(i>=0 and ps[i]*c<0)
  d8=add(*(scale(moments[h+l],vec[h]*vec[l])for h in range(5)for l in range(5)))
  h=add(scale(power(D,3),a0),scale(BN,a1))
  n=scale(d8,120)
  if i>=0:
   h=add(h,scale(mul(ys[i],power(D,2)),-c))
   n=add(n,scale(mul(P(ys[i]),power(D,2)),-c))
  if j>=0:
   bound=scale(Q(ys[j]),-6);direction=qs[j];bm=[(3*j+1,qs[j])]if j<6 else []
  else:
   j1,j2=test['b_pair'];assert 0<=j1<j2<8 and ps[j1]!=ps[j2]
   rank=[z+(z>=iv)+(z>=jv)for z in range(6)]+[iv,jv+1]
   direction=1 if ps[j1]*(rank[j1]-rank[j2])>0 else -1
   Y,Z=ys[j1],ys[j2]
   sums=lambda n:add(*(mul(power(Y,n-l),power(Z,l))for l in range(n+1)))
   bound=scale(add(sums(5),scale(mul(power(D,2),sums(3)),-15),scale(mul(BN,sums(2)),20),scale(mul(A4,add(Y,Z)),15)),-1)
   bm=[(3*z,ps[z])for z in(j1,j2)if z<6]
  n=add(n,mul(h,bound))
  hn='h5_'+str(index)+'_condition';nn='h5_'+str(index)+'_value'
  out[hn]=scale(h,direction);out[nn]=n
  mask=[]
  if i>=0 and i<6:mask.append((3*i,ps[i]))
  mask.extend(bm)
  masks[hn]=mask;masks[nn]=mask
 return out,masks

def full_extra_system(pat,tests):
 sig,iv,uq,up,jv,vq,vp=key(pat)
 one={(0,0,0):1};K={(1,0,0):1};S={(0,1,0):1};T={(0,0,1):1}
 D=add(power(K,2),scale(K,-1),one);AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
 CN=add(scale(power(K,2),2),scale(K,-2),scale(one,-1));DN=scale(add(AN,CN),-1)
 BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
 ys=[AN,scale(D,-1),CN,{},D,DN];ends=[scale(D,-3)]+ys+[scale(D,3)]
 U=add(ends[iv],mul(S,add(ends[iv+1],scale(ends[iv],-1))))
 V=add(ends[jv],mul(T,add(ends[jv+1],scale(ends[jv],-1))))
 A4=mul(mul(U,V),add(power(U,2),mul(U,V),power(V,2),scale(power(D,2),-6)))
 ys.extend([U,V]);ps=list(sig[::3])+[up,vp];qs=list(sig[1::3])+[uq,vq]
 def Q(Y):return add(power(Y,5),scale(mul(power(D,2),power(Y,3)),-10),scale(mul(BN,power(Y,2)),10),scale(mul(A4,Y),5))
 def P(Y):return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20),scale(mul(A4,power(Y,2)),15))
 # D^10 times the Newton moments at B=C=0.
 den=[power(D,j)for j in range(11)]
 coeff=[one,{},scale(power(D,2),-15),scale(BN,20),scale(A4,15),{},{}];moments=[scale(one,6)]
 for n in range(1,11):
  rhs=add(*(mul(coeff[j],moments[n-j])for j in range(1,min(n-1,6)+1)))
  if n<=6:rhs=add(rhs,scale(coeff[n],n))
  moments.append(scale(rhs,-1))
 moments=[mul(m,den[10-j])for j,m in enumerate(moments)]
 rank=[z+(z>=iv)+(z>=jv)for z in range(6)]+[iv,jv+1]
 def bound(meta):
  j=meta['q_node'];assert type(j)is int and -1<=j<8
  if j>=0:return scale(Q(ys[j]),-6),qs[j],[(3*j+1,qs[j])]if j<6 else []
  j1,j2=meta['b_pair'];assert 0<=j1<j2<8 and ps[j1]!=ps[j2]
  direction=1 if ps[j1]*(rank[j1]-rank[j2])>0 else -1
  Y,Z=ys[j1],ys[j2];sums=lambda n:add(*(mul(power(Y,n-l),power(Z,l))for l in range(n+1)))
  bb=scale(add(sums(5),scale(mul(power(D,2),sums(3)),-15),scale(mul(BN,sums(2)),20),scale(mul(A4,add(Y,Z)),15)),-1)
  return bb,direction,[(3*z,ps[z])for z in(j1,j2)if z<6]
 out={};masks={}
 for index,test in enumerate(tests):
  vec=test['vector'];i=test['p_node'];assert len(vec)==6 and all(type(z)is int for z in vec)and any(vec)
  assert type(i)is int and 0<=i<8
  v0,v1,v2,v3,v4,v5=vec;aa2=72000*v5*v5
  a0=-1200*(v0*v5+v1*v4+v2*v3)-25200*(v2*v5+v3*v4)-486000*v4*v5
  ab=19200*(2*v3*v5+v4*v4)+720000*v5*v5
  c0=-720*(2*v1*v5+2*v2*v4+v3*v3)-14400*(2*v3*v5+v4*v4)-270000*v5*v5
  a4=add(scale(den[4],a0),scale(mul(BN,D),ab),scale(A4,32400*v4*v5))
  c4=add(scale(den[4],c0),scale(mul(BN,D),43200*v4*v5),scale(A4,18000*v5*v5))
  d10=add(*(scale(moments[h+l],vec[h]*vec[l])for h in range(6)for l in range(6)))
  lo,sl,ml=bound(test['lower_bound']);hi,sh,mh=bound(test['upper_bound']);bd,direction,mb=bound(test)
  assert sl==1 and sh==-1
  astar=add(scale(mul(D,a4),120),scale(add(lo,hi),aa2))
  h=add(astar,scale(mul(c4,ys[i]),-120))
  n=add(scale(d10,14400),scale(mul(lo,hi),-aa2),scale(mul(c4,P(ys[i])),-120),mul(h,bd))
  pre='h6_'+str(index);out[pre+'_c']=scale(c4,ps[i]);out[pre+'_condition']=scale(h,direction);out[pre+'_value']=n
  mask=ml+mh+mb+([(3*i,ps[i])]if i<6 else [])
  for name in[pre+'_c',pre+'_condition',pre+'_value']:masks[name]=mask
 return out,masks

def raw_system(pat,tests):
 raw,e,masks=base_system(pat)
 assert tests and len({len(t['vector'])for t in tests})==1
 extra,em=(extra_system if len(tests[0]['vector'])==5 else full_extra_system)(pat,tests)
 raw.update(extra);masks.update(em)
 return raw,e,masks

def required_mask(cert):
 base=dict(cert);base['polynomials']=[p for p in cert['polynomials']if not p['name'].startswith(('h5_','h6_'))]
 mask=dict(base_mask(base));sig=key(cert['pattern'])[0]
 for row in cert['polynomials']:
  if not row['name'].startswith(('h5_','h6_')):continue
  i=int(row['name'].split('_')[1]);test=cert['linear_tests'][i]
  if 0<=test['p_node']<6:mask[3*test['p_node']]=sig[3*test['p_node']]
  if 0<=test['q_node']<6:mask[3*test['q_node']+1]=sig[3*test['q_node']+1]
  if test['q_node']==-1:
   for j in test['b_pair']:
    if j<6:mask[3*j]=sig[3*j]
  if len(test['vector'])==6:
   for bd in[test['lower_bound'],test['upper_bound']]:
    j=bd['q_node']
    if 0<=j<6:mask[3*j+1]=sig[3*j+1]
    if j==-1:
     for z in bd['b_pair']:
      if z<6:mask[3*z]=sig[3*z]
 return tuple(sorted(mask.items()))

def check_one(cert):
 pat=key(cert['pattern']);raw,e,masks=raw_system(pat,cert['linear_tests']);raw.update({'E+':e,'E-':scale(e,-1)});degrees=cert['degrees'];bs=[];shapes=[];names=[]
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
  if type(tree)is dict:
   if 'clause'in tree:return {tree['e']}|{z['i']for z in tree['claims']}
   return {tree['i'],tree['e']}
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
  if type(tree)is dict and 'clause'in tree:
   assert set(tree)=={'clause','claims','e'} and names[tree['e']]=='E+'
   j=tree['clause'];assert type(j)is int and 0<=j<len(cert['linear_tests'])
   prefix='h'+str(len(cert['linear_tests'][j]['vector']))+'_'+str(j)
   suffix=['_condition','_value']if len(cert['linear_tests'][j]['vector'])==5 else['_c','_condition','_value']
   assert len(tree['claims'])==len(suffix)
   assert [names[z['i']]for z in tree['claims']]==[prefix+t for t in suffix]
   b=elevate(arrays[tree['e']],shapes[tree['e']])
   for z in tree['claims']:
    assert set(z)=={'i','lambda'}
    q=F(*z['lambda']);a=elevate(arrays[z['i']],shapes[z['i']])
    assert max(x+q*y for x,y in zip(a,b))<=0
   stats['leaves']+=1;return
  if type(tree)is dict:
   assert not names[tree['i']].startswith(('h5_','h6_'))
   assert set(tree)=={'i','e','lambda'} and names[tree['e']]=='E+' and names[tree['i']]not in('E+','E-')
   q=F(*tree['lambda']);a=elevate(arrays[tree['i']],shapes[tree['i']]);b=elevate(arrays[tree['e']],shapes[tree['e']])
   assert max(x+q*y for x,y in zip(a,b))<=0
   stats['leaves']+=1;return
  if type(tree)is int:
   assert 0<=tree<len(bs)and tree in arrays and not names[tree].startswith(('h5_','h6_'))
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

def check_identities():
 """Check the vector formula symbolically, treating all five entries as variables."""
 from hermite_criterion import power as pp
 dim=9;zero=(0,)*dim;one={zero:1}
 def var(i):
  e=[0]*dim;e[i]=1;return {tuple(e):1}
 b,A,B,C=[var(i)for i in range(4)];v=[var(i)for i in range(4,9)]
 co=[one,{},scale(one,-15),scale(b,20),scale(A,60),scale(B,120),scale(C,120)]
 mom=[scale(one,6)]
 for n in range(1,9):
  rhs=add(*(mul(co[j],mom[n-j])for j in range(1,min(n-1,6)+1)))
  if n<=6:rhs=add(rhs,scale(co[n],n))
  mom.append(scale(rhs,-1))
 Q=add(*(mul(mul(v[i],v[j]),mom[i+j])for i in range(5)for j in range(5)))
 d={m:c for m,c in Q.items()if m[2]==m[3]==0}
 a=add(scale(add(mul(v[1],v[4]),mul(v[2],v[3])),-1200),scale(mul(v[3],v[4]),-25200),scale(mul(b,pp(v[4],2,dim)),19200))
 c=scale(add(scale(mul(v[2],v[4]),2),pp(v[3],2,dim),scale(pp(v[4],2,dim),20)),-720)
 assert Q==add(mul(a,B),mul(c,C),d)
 print('EXACTLY VERIFIED the formal H5 test-vector identity in nine variables.',flush=True)


def check_certificates():
 check_identities()
 check_full_identities()
 data=json.loads((ROOT/'data/hermite_fiber_certificates.json').read_text());seen=set();leaves=0
 for cert in data['certificates']:
  p=key(cert['pattern']);assert p not in seen;seen.add(p);leaves+=check_one(cert)['leaves']
 print('EXACTLY VERIFIED',len(seen),'Hermite-fibre certificates and',leaves,'three-dimensional Bernstein boxes.',flush=True)
 return seen

_certified=None
def reasons(word):
 global _certified
 if _certified is None:
  cs=json.loads((ROOT/'data/hermite_fiber_certificates.json').read_text())['certificates']
  _certified=[(key(c['pattern'])[1:],required_mask(c),i)for i,c in enumerate(cs)]
 return ['hermite_fiber_'+str(i)for p in patterns(word)for tail,mask,i in _certified
         if p[1:]==tail and all(p[0][j]==sg for j,sg in mask)]

def check_full_identities():
 from hermite_criterion import power as pp
 dim=10;one={(0,)*dim:1}
 def var(i):
  e=[0]*dim;e[i]=1;return {tuple(e):1}
 b,A,B,C=[var(i)for i in range(4)];v=[var(i)for i in range(4,10)];pw=lambda p,n:pp(p,n,dim)
 co=[one,{},scale(one,-15),scale(b,20),scale(A,60),scale(B,120),scale(C,120)];mom=[scale(one,6)]
 for n in range(1,11):
  rhs=add(*(mul(co[j],mom[n-j])for j in range(1,min(n-1,6)+1)))
  if n<=6:rhs=add(rhs,scale(co[n],n))
  mom.append(scale(rhs,-1))
 Q=add(*(mul(mul(v[i],v[j]),mom[i+j])for i in range(6)for j in range(6)))
 d={m:c for m,c in Q.items()if m[2]==m[3]==0}
 a=add(scale(add(mul(v[0],v[5]),mul(v[1],v[4]),mul(v[2],v[3])),-1200),scale(add(mul(v[2],v[5]),mul(v[3],v[4])),-25200),scale(mul(v[4],v[5]),-486000),scale(mul(b,add(scale(mul(v[3],v[5]),2),pw(v[4],2))),19200),scale(mul(b,pw(v[5],2)),720000),scale(mul(A,mul(v[4],v[5])),129600))
 c=add(scale(add(scale(mul(v[1],v[5]),2),scale(mul(v[2],v[4]),2),pw(v[3],2)),-720),scale(add(scale(mul(v[3],v[5]),2),pw(v[4],2)),-14400),scale(pw(v[5],2),-270000),scale(mul(b,mul(v[4],v[5])),43200),scale(mul(A,pw(v[5],2)),72000))
 assert Q==add(scale(mul(pw(v[5],2),pw(B,2)),72000),mul(a,B),mul(c,C),d)
 print('EXACTLY VERIFIED the formal H6 quadratic-in-B, affine-in-C identity in ten variables.',flush=True)

if __name__=='__main__':
 check_certificates()
