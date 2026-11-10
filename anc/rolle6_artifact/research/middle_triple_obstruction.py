"""A short exact obstruction for three consecutive zeros of P''.

The proof uses polynomial identities and elementary bounds, not subdivision.
"""
from fractions import Fraction as F
from two_parameter_obstruction import add,scale,mul,power
from triple_inflection_obstruction import ONE,U,R,identity_check as check_geometry
from inflection_obstruction import patterns

def check_identities():
    check_geometry()
    H=add(power(U,2),U,scale(mul(U,R),2),scale(R,-2),scale(ONE,-1))
    L=mul(add(scale(U,2),ONE),add(U,scale(ONE,2)))
    T=add(power(U,4),scale(power(U,3),-1),scale(power(U,2),-6),scale(U,-1),ONE)
    K={(3,1):4,(2,1):6,(1,1):-6,(0,1):-4,
       (4,0):3,(3,0):6,(2,0):-1,(1,0):-4,(0,0):-1}
    assert K==add(mul(L,H),T)
    t=add(ONE,U,scale(R,2))
    def primitive(x):
        return add(scale(power(x,6),F(1,30)),
                   scale(mul(add(ONE,R),power(x,5)),F(-1,10)),
                   scale(mul(H,power(x,4)),F(-1,12)),
                   scale(mul(mul(U,t),power(x,3)),F(1,6)))
    lhs=scale(add(primitive(scale(U,-1)),mul(U,primitive(ONE))),60)
    assert lhs==scale(mul(mul(U,add(U,ONE)),K),-1)
    K1={}
    for (i,j),c in K.items():K1[(i,0)]=K1.get((i,0),0)+c
    K1={ij:c for ij,c in K1.items()if c}
    assert K1=={(4,0):3,(3,0):10,(2,0):5,(1,0):-10,(0,0):-5}
    assert add(K,scale(K1,-1))==scale(mul(mul(add(R,scale(ONE,-1)),add(U,scale(ONE,-1))),L),2)
    half=F(1,2)
    small_bound=3*half**4+10*half**3+5*half**2-5
    assert small_bound==F(-37,16)
    dT={(i-1,j):c*i for (i,j),c in T.items()if i}
    positive=add(scale(mul(power(U,2),add(ONE,scale(U,-1))),4),
                 mul(U,add(scale(ONE,12),scale(U,-1))),ONE)
    assert scale(dT,-1)==positive
    assert sum(c*half**i for (i,j),c in T.items())==F(-17,16)
    # H(1,r)=1; H>0 for u>=1, r>1 because all terms are then nonnegative.
    assert H==add(power(U,2),U,scale(ONE,-1),scale(mul(R,add(U,scale(ONE,-1))),2))
    print('EXACTLY VERIFIED consecutive-triple identities; bounds -37/16 and -17/16; no Bernstein rectangles required.',flush=True)

def reasons(w):
    ps=patterns(w);sig=ps[0][0]
    if sig[5]==1 and sig[11]==-1 and [p[3]for p in ps[:3]]==[-1,1,-1]:
        return ['consecutive_inflections']
    return []

if __name__=='__main__':check_identities()
