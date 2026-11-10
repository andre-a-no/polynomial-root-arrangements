"""Analytic comparison of outer-inflection and two reference secants.

See PROOF_OUTER_SECANTS.md. This checker verifies algebra and word mapping;
it is not a formal proof assistant for the analytic inequalities.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
from left_minimum_inflection_obstruction import add,scale,mul
from hermite_criterion import power
from inflection_pair_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
CONDITIONS={'cubic_secant':[(0,1),(15,-1)],'unit_secant':[(3,1),(12,-1)],
            'mean_derivative':[(10,-1)],'middle_cubic_derivative':[(7,-1)]}
# Ascending coefficients in r, t. The denominators are positive.
NAD={(7,2):14,(6,2):-15,(5,1):-168,(4,2):32,(4,1):84,(3,0):288,
     (2,1):60,(2,0):-576,(1,2):-352,(1,1):2112,(1,0):-3168,
     (0,2):240,(0,1):-1440,(0,0):2160}
NPM={(6,2):7,(4,1):-84,(3,2):16,(3,1):-48,(2,0):144,(1,1):128,
     (1,0):-384,(0,2):-176,(0,1):1056,(0,0):-1584}
NC={(7,2):14,(6,2):30,(5,1):-168,(4,2):32,(4,1):-456,(3,0):288,
    (2,2):-216,(2,1):816,(2,0):72,(1,2):-352,(1,1):2112,(1,0):-3168,
    (0,2):-480,(0,1):2880,(0,0):-4320}
GC=[8,-12,-36,64,-31,4]
GC_NEG_BERNSTEIN=[F(16),F(40),F(472,5),F(992,5),F(1712,5),F(360)]
J=[48,64,-264,-40,231,20,-66,-8,3]
J_NEG_BERNSTEIN=[F(1344),F(2760),F(38432,7),F(73800,7),F(681672,35),
                 F(241200,7),F(405688,7),F(91440),F(131664)]

def check_identities():
 dim=3;one={(0,0,0):1};r={(1,0,0):1};d={(0,1,0):1};t={(0,0,1):1};ri={(-1,0,0):1}
 pw=lambda x,n:power(x,n,dim)
 const=lambda n:scale(one,n)
 sub_t=lambda p:{(i,j+2*k,0):v for (i,j,k),v in p.items()}
 def lift(p):return {(i,0,j):v for (i,j),v in p.items()}
 nad=lift(NAD);npm=lift(NPM);nc=lift(NC);T=add(pw(r,4),scale(pw(r,2),-2),const(4))
 w=mul(r,d);b=add(pw(d,3),scale(d,-3))
 p=add(const(-3),scale(ri,-6),mul(pw(d,2),add(scale(ri,2),scale(pw(r,2),F(1,2)))))
 q=add(pw(w,2),const(-6),scale(p,-1));A=scale(mul(p,q),F(1,4))
 assert scale(mul(w,add(p,scale(q,-1))),F(1,4))==b
 def divided(sum_,prod):
  seq=[one,sum_]
  for n in range(2,6):seq.append(add(mul(sum_,seq[-1]),scale(mul(prod,seq[-2]),-1)))
  return add(scale(seq[5],F(1,120)),scale(seq[3],F(-1,8)),scale(mul(b,seq[2]),F(1,6)),scale(mul(A,sum_),F(1,2)))
 uv=divided(w,p);ad=divided(d,add(pw(d,2),const(-3)))
 assert scale(mul(pw(r,2),add(uv,scale(ad,-1))),480)==mul(d,sub_t(nad))
 assert scale(mul(r,add(uv,scale(b,F(-1,6)))),240)==mul(d,sub_t(npm))
 qc=mul(d,add(scale(pw(d,4),F(9,20)),scale(pw(d,2),-1),scale(A,-1)))
 assert scale(mul(pw(r,2),add(uv,scale(qc,-1))),480)==mul(d,sub_t(nc))
 assert ad==scale(mul(d,add(const(-81),scale(A,20),scale(pw(d,2),23))),F(1,40))
 # Inner-root discriminant and product P'''(u)P'''(v).
 disc=add(pw(w,2),scale(q,-4))
 assert mul(r,disc)==mul(add(r,const(-2)),add(const(12),scale(mul(pw(d,2),add(pw(r,2),scale(r,2),const(4))),-1)))
 K=add(pw(w,2),scale(p,-1),const(-3));M=add(b,scale(mul(w,p),-1))
 prodS=add(mul(pw(K,2),p),mul(mul(K,M),w),pw(M,2))
 f1=add(mul(add(pw(r,2),const(2)),pw(d,2)),const(-6))
 f2=add(mul(add(pw(r,2),scale(r,-2),const(4)),pw(d,2)),const(-12))
 f3=add(mul(T,pw(d,2)),const(-12))
 assert scale(mul(pw(r,3),prodS),8)==mul(mul(mul(add(r,const(2)),f1),f2),f3)
 # H3, written using independent r,t.
 hnum=add(mul(pw(r,6),pw(t,2)),scale(mul(pw(r,4),t),-12),scale(mul(pw(r,2),pw(t,3)),8),scale(mul(pw(r,2),pw(t,2)),-48),scale(mul(pw(r,2),t),72),scale(pw(r,2),16),scale(pw(t,2),-16),scale(t,96),const(-144))
 H=add(const(5),scale(A,-4),scale(pw(b,2),-2))
 assert scale(mul(pw(r,2),H),4)==scale(sub_t(hnum),-1)
 def endpoint(n,degree):
  return add(*(scale(mul({(i,0,0):1},pw(T,degree-j)),v*12**j)for (i,_,j),v in n.items()))
 # -H(r,12/T)*r^2*T^3: factored expression.
 eh=scale(endpoint(hnum,3),F(1,4))
 rhs=scale(mul(mul(mul(mul(pw(r,2),pw(add(r,const(-2)),2)),pw(add(r,const(-1)),2)),pw(add(r,const(1)),2)),mul(pw(add(r,const(2)),2),add(pw(r,4),scale(pw(r,2),-14),const(4)))),4)
 assert eh==rhs
 # Quadratic coefficients and endpoint identities for both secant differences.
 coeff=lambda n,k:{(i,0,0):v for (i,_,j),v in n.items()if j==k}
 assert coeff(nad,2)==mul(add(pw(r,3),const(-4)),add(scale(pw(r,4),14),scale(pw(r,3),-15),scale(r,88),const(-60)))
 assert coeff(nad,0)==scale(add(scale(pw(r,3),2),scale(pw(r,2),-4),scale(r,-22),const(15)),144)
 G=add(pw(r,5),scale(pw(r,4),2),scale(pw(r,3),-20),scale(pw(r,2),-69),scale(r,-72),const(-22))
 assert endpoint(nad,2)==scale(mul(mul(pw(r,2),pw(add(r,const(-1)),4)),G),288)
 assert coeff(npm,2)==mul(add(pw(r,3),const(-4)),add(scale(pw(r,3),7),const(44)))
 assert coeff(npm,0)==scale(add(scale(pw(r,2),3),scale(r,-8),const(-33)),48)
 assert endpoint(npm,2)==scale(mul(pw(r,2),{(i,0,0):v for i,v in enumerate(J)if v}),48)
 assert coeff(nc,2)==scale(mul(pw(add(r,const(2)),2),add(scale(pw(r,5),7),scale(pw(r,4),-13),scale(pw(r,3),24),scale(pw(r,2),-28),scale(r,16),const(-60))),2)
 assert coeff(nc,0)==scale(mul(pw(add(r,const(2)),2),add(scale(r,4),const(-15))),72)
 assert endpoint(nc,2)==scale(mul(mul(pw(r,2),pw(add(r,const(2)),4)),{(i,0,0):v for i,v in enumerate(GC)if v}),72)
 # Exact Bernstein coefficients of -J(2+2z), degree 8.
 qcoef=[-sum(F(J[k]*comb(k,j)*2**k)for k in range(j,9))for j in range(9)]
 bs=[sum(qcoef[j]*F(comb(i,j),comb(8,j))for j in range(i+1))for i in range(9)]
 assert bs==J_NEG_BERNSTEIN and min(bs)>0
 gc_power=[-sum(F(GC[k]*comb(k,j)*2**k)for k in range(j,6))for j in range(6)]
 gc_bs=[sum(gc_power[j]*F(comb(i,j),comb(5,j))for j in range(i+1))for i in range(6)]
 assert gc_bs==GC_NEG_BERNSTEIN and min(gc_bs)>0
 assert F(675,16)-F(1350,193)-6-18-F(128,75)==F(2196961,231600)>0
 assert F(15,4)**4-14*F(15,4)**2+4==F(1249,256)>0
 # H_t>0 on r>=4, 0<=t<=12/T: the lower bound used in the proof.
 assert F(48)-18-F(48,7)-6-F(3,2)==F(219,14)>0
 # Independent 3x3 Newton moment determinant, variables A,b in two slots.
 aa={(1,0,0):1};bb={(0,1,0):1}
 det=add(scale(add(const(13500),scale(aa,-7200),scale(pw(bb,2),-3600)),6),const(-27000))
 assert det==scale(add(const(5),scale(aa,-4),scale(pw(bb,2),-2)),10800)
 print('EXACTLY VERIFIED outer-secant identities, H3 determinant, and fifteen positive Bernstein coefficients.',flush=True)
 return {'bernstein_coefficients':15,'new_rules':4}

def reasons(word):
 out=[]
 for sig,iv,uq,up,jv,vq,vp in patterns(word):
  if (iv,jv,up,vp)!=(0,6,-1,1):continue
  for name,mask in CONDITIONS.items():
   if all(sig[i]==sg for i,sg in mask):out.append('outer_secant_'+name)
 return sorted(set(out))

if __name__=='__main__':check_identities()
