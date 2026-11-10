"""Exact weighted-average inequality at outer inflections of P-prime."""
from fractions import Fraction as F
from pathlib import Path
import json
from two_parameter_obstruction import add,scale,mul,power
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
ONE={(0,0):1};M={(1,0):1};P={(0,1):1}

def check_identities():
 e=scale(mul(M,P),2)
 t=mul(M,add(ONE,scale(P,-3)))
 w=add(ONE,scale(P,-2),scale(mul(M,t),-2))
 k2=add(power(e,2),scale(w,-1))
 assert k2==mul(add(scale(P,2),scale(ONE,-1)),add(ONE,scale(mul(power(M,2),add(ONE,scale(P,-1))),-2)))
 Rc=mul(add(power(M,2),scale(ONE,-1)),add(power(add(e,scale(M,-1)),2),scale(k2,-1)))
 assert Rc==mul(add(scale(P,2),scale(ONE,-1)),power(add(ONE,scale(power(M,2),-1)),2))
 S2=add(power(t,2),P);S3=add(power(t,3),scale(mul(P,t),2))
 intercept=add(mul(P,add(scale(S3,F(1,5)),scale(mul(e,S2),F(1,2)),
                           scale(mul(add(w,scale(ONE,-1)),t),F(1,3)),scale(e,-1))),scale(e,F(7,30)))
 positive=add(mul(mul(power(M,2),power(add(scale(P,3),scale(ONE,-1)),2)),add(scale(ONE,7),scale(P,-6))),
              mul(add(ONE,scale(P,-1)),add(scale(P,27),scale(ONE,-7))))
 assert intercept==scale(mul(mul(M,P),positive),F(-1,15))
 assert F(1,2)*27-7>0 and 7-6*F(1)>0
 print('EXACTLY VERIFIED outer-inflection weighted-average identity; strictly positive remainder for 0<m<1, 1/2<p<1.',flush=True)

def reasons(w):
 ps=patterns(w);sig=ps[0][0]
 if sig[1]==1 and sig[16]==1 and ps[0][3]==1 and ps[3][3]==-1:
  return ['outer_inflection_average']
 return []

if __name__=='__main__':
 check_identities()
 known=json.loads((ROOT/'data/certificates_n6.json').read_text());unknown=json.loads((ROOT/'data/unresolved_n6.json').read_text())
 assert not any(reasons(w)for w in known)
 print('newly matched:',[w for w in unknown if reasons(w)])
