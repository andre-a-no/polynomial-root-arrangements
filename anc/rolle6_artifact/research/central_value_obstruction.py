"""Exact affine identities at the roots of p^(4) and their midpoint.

Normalize the roots of p^(4) to -1,+1 and p''' to x^3-3x+b.
Reflect if necessary to have b>0. At -1,0,+1 the values of p,p',p''
are affine in b,A,B,C. A signed sum of strictly positive values cannot
be a nonpositive constant. See PROOF_CENTRAL_VALUES_AND_OVERLAP.md.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
# Each row is 60 times the coefficients of b,A,B,C and the constant.
# Order: p(-1),p'(-1),p''(-1),p(0),p'(0),p''(0),p(1),p'(1),p''(1),b.
FORMS60=[
    [-10,30,-60,60,-7], [30,-60,60,0,27], [-60,60,0,0,-75],
    [0,0,0,60,0], [0,0,60,0,0], [0,60,0,0,0],
    [10,30,60,60,-7], [30,60,60,0,-27], [60,60,0,0,-75],
    [60,0,0,0,0],
]
CERTIFICATES=json.loads((ROOT/'data/central_value_identities.json').read_text())


def canonical(w):
    t=[j for j,c in enumerate(w)if c=='3']
    return w if t[1]>w.index('5')else w[::-1]


def signature(w):
    w=canonical(w)
    assert len(w)==21 and all(w.count(str(i))==6-i for i in range(6))
    pos=[[j for j,c in enumerate(w)if c==str(i)]for i in range(6)]
    markers=sorted(pos[4]+pos[5])
    return tuple((-1)**sum(j>x for j in pos[d])for x in markers for d in range(3))


def check_identities():
    from fractions import Fraction as F
    from math import factorial
    # Independently differentiate the normalized sextic and evaluate it.
    coefficients={6:[0,0,0,0,F(1,120)],4:[0,0,0,0,F(-1,8)],
        3:[F(1,6),0,0,0,0],2:[0,F(1,2),0,0,0],
        1:[0,0,1,0,0],0:[0,0,0,1,0]}
    generated=[]
    for x in (-1,0,1):
        for d in range(3):
            generated.append([60*sum(v[j]*F(factorial(k),factorial(k-d))*x**(k-d)
                for k,v in coefficients.items()if k>=d)for j in range(5)])
    assert generated==FORMS60[:9]
    assert len({c['id']for c in CERTIFICATES})==len(CERTIFICATES)
    assert len({tuple(c['weights'])for c in CERTIFICATES})==len(CERTIFICATES)
    for c in CERTIFICATES:
        w=c['weights']
        assert len(w)==10 and all(isinstance(a,int)for a in w)
        assert w[-1]>=0 and any(w)
        total=[sum(a*f[j]for a,f in zip(w,FORMS60))for j in range(5)]
        assert total[:4]==[0]*4
        assert total[-1]==c['constant_times_60']<=0
    print('EXACTLY VERIFIED',len(CERTIFICATES),'distinct central-value identities',flush=True)


def reasons(w):
    signs=signature(w)
    return [c['id']for c in CERTIFICATES
        if all(a==0 or a*s>0 for a,s in zip(c['weights'][:9],signs))]


if __name__=='__main__':
    check_identities()
