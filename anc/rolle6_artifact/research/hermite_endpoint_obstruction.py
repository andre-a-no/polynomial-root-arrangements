"""A marker-sign obstruction proved by convexity and two univariate signs.

For c<0, six simple real roots forbid P(-1)>0, P'(a)<0<P'(d).
All algebraic identities and the final two Bernstein sign tests are exact.
"""
import json,math
from fractions import Fraction as F
from pathlib import Path
from left_minimum_inflection_obstruction import add,scale,mul
from hermite_criterion import power
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]

def check_identities():
 # Independently reconstruct Q_v from Newton sums.
 dim=4;one={(0,)*dim:1}
 def var(i):
  m=[0]*dim;m[i]=1;return {tuple(m):1}
 b,A,B,C=[var(i)for i in range(4)];pw=lambda p,n:power(p,n,dim)
 co=[one,{},scale(one,-15),scale(b,20),scale(A,60),scale(B,120),scale(C,120)];mom=[scale(one,6)]
 for n in range(1,9):
  rhs=add(*(mul(co[j],mom[n-j])for j in range(1,min(n-1,6)+1)))
  if n<=6:rhs=add(rhs,scale(co[n],n))
  mom.append(scale(rhs,-1))
 v=[-7,0,8,0,-1];Q=add(*(scale(mom[i+j],v[i]*v[j])for i in range(5)for j in range(5)))
 dd=add(scale(pw(A,2),14400),scale(A,-40320),scale(pw(b,2),28800),scale(one,25284))
 assert Q==add(scale(mul(b,B),19200),scale(C,-2880),dd)
 # Work modulo d^2+c*d+c^2-3=0, with a=-c-d and b=3c-c^3.
 dim=3;one={(0,0,0):1};c={(1,0,0):1};d={(0,1,0):1};A={(0,0,1):1};pw=lambda p,n:power(p,n,dim)
 a=scale(add(c,d),-1);b=add(scale(c,3),scale(pw(c,3),-1));rel=add(scale(one,3),scale(pw(c,2),-1),scale(mul(c,d),-1))
 def reduce(p):
  p=dict(p)
  while any(m[1]>=2 for m in p):
   m=max((m for m in p if m[1]>=2),key=lambda m:m[1]);q=p.pop(m);e=(m[0],m[1]-2,m[2]);p=add(p,scale(mul({e:1},rel),q))
  return p
 T=scale(add(scale(one,21),scale(pw(c,2),27),scale(pw(c,4),-9)),F(1,20))
 U=scale(add(scale(one,5),scale(pw(b,2),-2)),F(1,4))
 q0=lambda x,AA:add(scale(pw(x,5),F(1,20)),scale(pw(x,3),F(-1,2)),scale(mul(b,pw(x,2)),F(1,2)),mul(AA,x))
 star=scale(mul(b,add(scale(one,34),scale(pw(c,2),-9))),F(1,20))
 assert reduce(add(q0(d,A),scale(q0(a,A),-1)))==reduce(mul(add(d,scale(a,-1)),add(A,scale(T,-1))))
 assert reduce(q0(d,T))==star
 assert scale(add(U,scale(T,-1)),20)==add(scale(one,4),scale(pw(c,2),-117),scale(pw(c,4),69),scale(pw(c,6),-10))
 # N(A)>0 is forced by Q_v>0, P(-1)>0 and P'(d)>0.
 g=lambda AA:add(scale(pw(AA,2),14400),scale(AA,-38880),scale(pw(b,2),28800),scale(b,-480),scale(one,24948))
 coef=add(scale(one,2880),scale(b,-19200))
 NT=add(g(T),mul(coef,star))
 R=add(scale(pw(c,5),963),scale(pw(c,3),-3209),scale(pw(c,2),108),scale(c,1284),scale(one,-368))
 assert NT==scale(mul(mul(c,add(pw(c,2),scale(one,-3))),R),12)
 M=scale(one,F(11,10));upper=add(g(M),mul(coef,add(star,scale(add(M,scale(T,-1)),2))))
 expected=scale(add(scale(pw(c,8),720),scale(pw(c,7),1440),scale(pw(c,6),-4640),scale(pw(c,5),-8532),scale(pw(c,4),8616),scale(pw(c,3),12428),scale(pw(c,2),-3528),scale(c,624),scale(one,-9)),12)
 assert upper==expected
 # R0(t)=4-117t+69t^2-10t^3 decreases on [0,1].
 t=F(4,49);assert 4-117*t+69*t*t-10*t*t*t<0
 # Its derivative equals -9-(1-t)(108-30t)<0 on that interval.
 assert [-117,138,-30]==[-9-108,108+30,-30]
 # Hermite F(-1) at M equals b(464b+728)-76/5<0 for -1<b<0.
 assert -1600*F(11,10)**2+1728*F(11,10)+20==F(-76,5)
 assert 480*F(11,10)+200==728
 data=json.loads((ROOT/'data/hermite_endpoint_certificate.json').read_text())
 assert data['vector']==v and data['marker_conditions']==[[1,-1],[3,1],[16,1]]
 for name,p in [('minus_R',scale(R,-1)),('minus_upper',scale(upper,-1))]:
  assert all(j==h==0 for i,j,h in p)
  n=max(i for i,j,h in p);co=[F(p.get((i,0,0),0))*F(-2,7)**i for i in range(n+1)]
  bs=[sum((co[j]*F(math.comb(i,j),math.comb(n,j))for j in range(i+1)),F(0))for i in range(n+1)]
  assert all(q>0 for q in bs)
  row=data[name];assert row['power_coefficients']==[[q.numerator,q.denominator]for q in co]
  assert row['bernstein_coefficients']==[[q.numerator,q.denominator]for q in bs]
 print('EXACTLY VERIFIED Hermite endpoint obstruction: Newton identities, convexity endpoints and 15 positive univariate Bernstein coefficients.',flush=True)

def reasons(word):
 sig=patterns(word)[0][0]
 return ['hermite_endpoint']if sig[1]==-1 and sig[3]==1 and sig[16]==1 else[]

if __name__=='__main__':check_identities()
