"""An analytic root-order obstruction for sextics, with exact identity checks.

If s1<s2<s3<s4 are the roots of p'' and t1<t2<t3 the roots of p''',
then (p(s2)-p(s1))/(s2-s1) < p'(t1) and
(p(s4)-p(s3))/(s4-s3) > p'(t3), assuming positive leading coefficient.
The proof is in PROOF_FIRST_FAMILY.md. Only the standard library is used.
"""
from fractions import Fraction

def add(*polys):
    out={}
    for p in polys:
        for m,c in p.items():out[m]=out.get(m,0)+c
    return {m:c for m,c in out.items()if c}
def scale(p,k):return {m:k*c for m,c in p.items()if k*c}
def mul(p,q):
    out={}
    for u,a in p.items():
        for v,b in q.items():
            m=tuple(i+j for i,j in zip(u,v));out[m]=out.get(m,0)+a*b
    return {m:c for m,c in out.items()if c}
def power(p,n):
    out={(0,0,0):1}
    for _ in range(n):out=mul(out,p)
    return out
def derivative(p):return {(m[0]-1,m[1],m[2]):m[0]*c for m,c in p.items()if m[0]}
def evaluate_univariate(p,x):
    assert all(m[1:]==(0,0)for m in p)
    return sum(c*x**m[0]for m,c in p.items())

def check_identities():
    one={(0,0,0):1};u={(1,0,0):1};S={(0,1,0):1};T={(0,0,1):1}
    A={(4,0,0):15,(3,0,0):-30,(2,0,0):5,(1,0,0):10,(0,0,0):-2}
    B={(5,0,0):-28,(4,0,0):85,(3,0,0):-70,(2,0,0):-5,(1,0,0):20,(0,0,0):-3}
    E60=add(scale(one,-1),scale(S,2),scale(T,-5),scale(power(u,5),-12),
        scale(mul(add(one,S),power(u,4)),15),scale(mul(add(S,T),power(u,3)),-20),scale(mul(T,power(u,2)),30))
    # Independently integrate R(x)=x(x-1)(x*x-S*x+T) term by term.
    rcoeff=[{},scale(T,-1),add(S,T),scale(add(one,S),-1),one]
    integrated_difference={}
    for k,coeff in enumerate(rcoeff):
        integrated_difference=add(integrated_difference,
            scale(coeff,Fraction(1,(k+1)*(k+2))),
            scale(mul(coeff,power(u,k+1)),Fraction(-1,k+1)))
    assert scale(integrated_difference,60)==E60
    Rprime=add(scale(power(u,3),4),scale(mul(add(one,S),power(u,2)),-3),scale(mul(add(S,T),u),2),scale(T,-1))
    M=add(scale(power(u,2),-10),scale(u,10),scale(one,5))
    assert not add(E60,mul(add(S,scale(one,-2)),A),B,scale(mul(M,Rprime),-1))
    assert derivative(A)==scale(mul(add(scale(u,2),scale(one,-1)),add(scale(power(u,2),3),scale(u,-3),scale(one,-1))),10)
    v=add(scale(u,4),scale(one,-1))
    bracket=add(scale(power(v,3),7),scale(power(v,2),-30),scale(v,3),scale(one,48))
    assert not add(scale(derivative(derivative(B)),4),scale(bracket,5))
    assert evaluate_univariate(A,Fraction(1,4))==Fraction(103,256)
    assert evaluate_univariate(B,Fraction(1,4))==Fraction(115,128)
    assert evaluate_univariate(B,Fraction(1,2))==Fraction(23,16)
    print('EXACTLY VERIFIED analytic mean-value identity, derivative identities and endpoint bounds',flush=True)

def reasons(w):
    assert len(w)==21 and all(w.count(str(i))==6-i for i in range(6))
    pos=[[j for j,c in enumerate(w)if c==str(i)]for i in range(6)]
    sign=lambda row,x:(-1)**sum(j>x for j in pos[row])
    out=[]
    if sign(0,pos[2][0])<0 and sign(0,pos[2][1])>0 and sign(1,pos[3][0])<0:out.append('first_interval')
    if sign(0,pos[2][2])>0 and sign(0,pos[2][3])<0 and sign(1,pos[3][2])>0:out.append('last_interval')
    return out

if __name__=='__main__':check_identities()
