"""Exact Hermite tests and the identities used by the pair obstruction."""
from fractions import Fraction as F
from itertools import permutations
from left_minimum_inflection_obstruction import add,scale,mul

def power(p,n,dim):
 q={(0,)*dim:1}
 for _ in range(n):q=mul(q,p)
 return q

def det_poly(matrix,dim):
 n=len(matrix);out={}
 for perm in permutations(range(n)):
  sg=(-1)**sum(perm[i]>perm[j]for i in range(n)for j in range(i+1,n));p={(0,)*dim:sg}
  for i,j in enumerate(perm):p=mul(p,matrix[i][j])
  out=add(out,p)
 return out

def check_identities():
 one={(0,0,0,0):1};b={(1,0,0,0):1};A={(0,1,0,0):1};B={(0,0,1,0):1};C={(0,0,0,1):1}
 pw=lambda p,n:power(p,n,4)
 coeff=[one,{},scale(one,-15),scale(b,20),scale(A,60),scale(B,120),scale(C,120)]
 moments=[scale(one,6)]
 for k in range(1,7):
  moments.append(scale(add(*(mul(coeff[j],moments[k-j])for j in range(1,k)),scale(coeff[k],k)),-1))
 d3=det_poly([[moments[i+j]for j in range(3)]for i in range(3)],4)
 G=add(scale(one,5),scale(A,-4),scale(pw(b,2),-2))
 assert d3==scale(G,10800)
 d4=det_poly([[moments[i+j]for j in range(4)]for i in range(4)],4)
 K0=add(scale(pw(A,3),128),scale(pw(A,2),-280),scale(mul(A,pw(b,2)),-180),scale(A,150),scale(pw(b,4),-20),scale(pw(b,2),25))
 formula=add(K0,scale(pw(B,2),-100),scale(mul(mul(add(scale(A,160),scale(one,100)),b),B),-1),scale(mul(G,C),-12))
 assert d4==scale(formula,648000)
 # Eliminate C at a required positive value P(x), then maximize over B.
 one={(0,0,0):1};b={(1,0,0):1};A={(0,1,0):1};x={(0,0,1):1};pw=lambda p,n:power(p,n,3)
 G=add(scale(one,5),scale(A,-4),scale(pw(b,2),-2))
 K0=add(scale(pw(A,3),128),scale(pw(A,2),-280),scale(mul(A,pw(b,2)),-180),scale(A,150),scale(pw(b,4),-20),scale(pw(b,2),25))
 g0=add(scale(pw(x,6),F(1,120)),scale(pw(x,4),F(-1,8)),scale(mul(b,pw(x,3)),F(1,6)),scale(mul(A,pw(x,2)),F(1,2)))
 L=add(scale(mul(G,x),12),scale(mul(add(scale(A,160),scale(one,100)),b),-1))
 K=add(K0,scale(mul(G,g0),12))
 h=add(scale(pw(A,2),-1600),scale(mul(mul(A,b),x),-480),scale(mul(A,pw(x,2)),228),scale(A,1500),
       scale(mul(pw(b,2),pw(x,2)),-36),scale(pw(b,2),500),scale(mul(b,pw(x,3)),100),scale(mul(b,x),-300),
       scale(pw(x,6),5),scale(pw(x,4),-75),scale(pw(x,2),90))
 assert add(pw(L,2),scale(K,400))==scale(mul(G,h),8)
 print('EXACTLY VERIFIED Hermite Delta3, Delta4 and the elimination identity L^2+400K=8GF.',flush=True)

def moments(coefficients,maximum):
 """Monic coefficients in descending degree order, exact rational input."""
 a=list(map(F,coefficients));a=[v/a[0]for v in a];n=len(a)-1;s=[F(n)]
 for k in range(1,maximum+1):
  s.append(-sum((a[j]*s[k-j]for j in range(1,min(k-1,n)+1)),F(0))-(k*a[k]if k<=n else 0))
 return s

def principal_minors(coefficients):
 n=len(coefficients)-1;s=moments(coefficients,2*n-2);out=[]
 for size in range(1,n+1):
  m=[[s[i+j]for j in range(size)]for i in range(size)];sign=1
  for i in range(size):
   pivot=next((j for j in range(i,size)if m[j][i]),None)
   if pivot is None:sign=0;break
   if pivot!=i:m[i],m[pivot]=m[pivot],m[i];sign=-sign
   for j in range(i+1,size):
    q=m[j][i]/m[i][i]
    for h in range(i,size):m[j][h]-=q*m[i][h]
  v=F(sign)
  for i in range(size):v*=m[i][i]
  out.append(v)
 return out

def check_audits():
 import json
 from pathlib import Path
 from verify import polynomial
 root=Path(__file__).resolve().parents[1];folder=root/'research/work/inflection_pairs'
 data=json.loads((folder/'hermite_audit.json').read_text());known=json.loads((root/'data/certificates_n6.json').read_text())
 assert sorted(known)==data['known_words'] and len(known)==data['known_hyperbolic_words_accepted']
 for item in known.values():
  coeff=polynomial(item['roots'])if 'roots'in item else item['coefficients']
  assert all(v>0 for v in principal_minors(coeff))
 example=json.loads((root/'data/single_stationary_limit_example.json').read_text())
 bad=principal_minors(example['constructed_coefficients'])
 assert [[v.numerator,v.denominator]for v in bad]==data['nonhyperbolic_stationary_example_minors']
 assert [(v>0)-(v<0)for v in bad]==data['nonhyperbolic_stationary_example_minor_signs']
 assert not all(v>0 for v in bad)
 examples=json.loads((root/'data/birkhoff_eight_point_audit.json').read_text())['new_witnesses']
 rows=[]
 for name,item in examples.items():
  values=principal_minors(item['coefficients']);assert not all(v>0 for v in values)
  rows.append({'configuration':name,'minor_signs':[(v>0)-(v<0)for v in values],'real_root_count':item['exact_real_root_counts'][0]})
 assert rows==json.loads((folder/'hermite_interpolation_audit.json').read_text())
 print('EXACTLY VERIFIED full Hermite audits:',len(known),'hyperbolic witnesses accepted; 53 nonhyperbolic witnesses rejected.',flush=True)

if __name__=='__main__':
 check_identities()
 check_audits()
