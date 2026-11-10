"""Two analytic sign obstructions from Laguerre and an integral identity.

The analytic argument is in PROOF_LAGUERRE_TRAPEZOID.md. Only the exact
polynomial identities and mapping to words are checked by this module.
"""
from pathlib import Path
from fractions import Fraction as F
from left_minimum_inflection_obstruction import add,scale,mul
from hermite_criterion import power
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
CONDITIONS={'outer_cubic':[(0,-1),(1,1),(15,-1),(16,1)],
            'fourth_derivative':[(3,-1),(4,1),(5,-1),(12,-1),(13,1)]}

def check_identities():
 # Variables a,c,d,A,B,C. Impose a+c+d=0 and b=-a*c*d.
 dim=6;one={(0,)*dim:1}
 def var(i):
  e=[0]*dim;e[i]=1;return {tuple(e):1}
 a,c,d,A,B,C=[var(i)for i in range(6)];pw=lambda p,n:power(p,n,dim)
 # Integrate the general cubic (x-a)(x-c)(x-d) three times.
 ss=add(a,c,d);pairs=add(mul(a,c),mul(a,d),mul(c,d));prod=mul(mul(a,c),d)
 def P(x):return add(scale(pw(x,6),F(1,120)),scale(mul(ss,pw(x,5)),F(-1,60)),scale(mul(pairs,pw(x,4)),F(1,24)),scale(mul(prod,pw(x,3)),F(-1,6)),scale(mul(A,pw(x,2)),F(1,2)),mul(B,x),C)
 def Q(x):return add(scale(pw(x,5),F(1,20)),scale(mul(ss,pw(x,4)),F(-1,12)),scale(mul(pairs,pw(x,3)),F(1,6)),scale(mul(prod,pw(x,2)),F(-1,2)),mul(A,x),B)
 h=add(d,scale(a,-1));mid_minus_c=add(scale(add(a,d),F(1,2)),scale(c,-1))
 rhs=add(scale(mul(h,add(Q(a),Q(d))),F(1,2)),scale(mul(pw(h,5),mid_minus_c),F(1,60)))
 assert add(P(d),scale(P(a),-1))==rhs
 # If a+c+d=0, midpoint-c=-3c/2 and the correction is -c*h^5/40.
 assert add(mid_minus_c,scale(c,F(3,2)))==scale(ss,F(1,2))
 # Direct normalized endpoint identity at -1,+1.
 dim=4;one={(0,)*dim:1};b={(1,0,0,0):1};A={(0,1,0,0):1};B={(0,0,1,0):1};C={(0,0,0,1):1}
 def P(x):return add(scale(one,F(x**6,120)-F(x**4,8)),scale(b,F(x**3,6)),scale(A,F(x*x,2)),scale(B,x),C)
 def Q(x):return add(scale(one,F(x**5,20)-F(x**3,2)),scale(b,F(x*x,2)),scale(A,x),B)
 assert add(P(1),scale(P(-1),-1))==add(Q(-1),Q(1),scale(b,F(-2,3)))
 # Reflection is handled by the existing canonical marker projection.
 print('EXACTLY VERIFIED both polynomial integral identities for the Laguerre-trapezoid obstructions.',flush=True)

def reasons(word):
 sig=patterns(word)[0][0]
 return ['laguerre_trapezoid_'+name for name,mask in CONDITIONS.items()if all(sig[i]==sg for i,sg in mask)]

if __name__=='__main__':check_identities()
