"""Polynomial Farkas certificates covering the full parameter interval (0,1).

Normalize p^(4) roots to -1,1, reflect so its middle p''' root v lies
in (0,1), and substitute b=3v-v^3. The constraints are affine in A,B,C.
Polynomial multipliers and their signs are checked with exact Bernstein
coefficients on rational subintervals. No external packages are needed.
"""
import json
from fractions import Fraction as F
from math import comb,factorial
from pathlib import Path
from central_value_obstruction import canonical

ROOT=Path(__file__).resolve().parents[1]
DATA=json.loads((ROOT/'data/parametric_middle_certificates.json').read_text())
# Coefficients of A,B,C,constant, each a polynomial in v (low powers first).
# Values are scaled by 120; rows are p,p',p'' at -1,0,v,1, in that order.
FORMS=[
    [[60],[-120],[120],[-14,-60,0,20]],
    [[-120],[120],[],[54,180,0,-60]],
    [[120],[],[],[-150,-360,0,120]],
    [[],[],[120],[]],[[],[120],[],[]],[[120],[],[],[]],
    [[0,0,60],[0,120],[120],[0,0,0,0,45,0,-19]],
    [[0,120],[120],[],[0,0,0,120,0,-54]],
    [[120],[],[],[0,0,180,0,-90]],
    [[60],[120],[120],[-14,60,0,-20]],
    [[120],[120],[],[-54,180,0,-60]],
    [[120],[],[],[-150,360,0,-120]],
]


def trim(p):
    p=list(p)
    while p and not p[-1]:p.pop()
    return p


def add(*polys):
    out=[0]*max((len(p)for p in polys),default=0)
    for p in polys:
        for i,c in enumerate(p):out[i]+=c
    return trim(out)


def scale(p,c):return trim([c*x for x in p])


def mul(p,q):
    if not p or not q:return []
    out=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):out[i+j]+=a*b
    return trim(out)


def power(p,n):
    out=[1]
    for _ in range(n):out=mul(out,p)
    return out


def value(p,x):
    out=0
    for c in reversed(p):out=out*x+c
    return out


def bernstein(p,left,right):
    p=trim(p)
    if not p:return [F(0)]
    n=len(p)-1
    # Coefficients after v=left+(right-left)*t, 0<=t<=1.
    a=[sum(p[i]*comb(i,k)*left**(i-k)*(right-left)**k
        for i in range(k,n+1))for k in range(n+1)]
    return [sum(a[i]*F(comb(k,i),comb(n,i))for i in range(k+1))for k in range(n+1)]


def signature(w):
    w=canonical(w)
    assert len(w)==21 and all(w.count(str(i))==6-i for i in range(6))
    pos=[[j for j,c in enumerate(w)if c==str(i)]for i in range(6)]
    markers=sorted(pos[4]+pos[5]+[pos[3][1]])
    return tuple((-1)**sum(j>x for j in pos[d])for x in markers for d in range(3))


MASKS={tuple(p['signature']):p['id']for p in DATA['patterns']}


def reasons(w):
    signs=signature(w)
    return [name for mask,name in MASKS.items()
        if all(a==0 or a==b for a,b in zip(mask,signs))]


def check_certificates():
    # Independently differentiate the normalized sextic at polynomial nodes.
    xcoeff={6:[[],[],[],[F(1,120)]],4:[[],[],[],[F(-1,8)]],
        3:[[],[],[],[0,F(1,2),0,F(-1,6)]],
        2:[[F(1,2)],[],[],[]],1:[[],[1],[],[]],0:[[],[],[1],[]]}
    generated=[]
    for x in ([-1],[],[0,1],[1]):
        for d in range(3):
            row=[]
            for j in range(4):
                row.append(add(*(scale(mul(cs[j],power(x,k-d)),
                    120*factorial(k)//factorial(k-d))for k,cs in xcoeff.items()if k>=d)))
            generated.append(row)
    assert generated==FORMS
    identities={c['id']:c for c in DATA['identities']}
    assert len(identities)==len(DATA['identities'])
    assert len(MASKS)==len(DATA['patterns'])
    # Zero means an unconstrained sign. No retained mask implies another.
    for a in MASKS:
        for b in MASKS:
            if a!=b:assert not all(x==0 or x==y for x,y in zip(a,b))
    assert len({p['id']for p in DATA['patterns']})==len(DATA['patterns'])
    for cert in identities.values():
        terms=cert['terms']
        assert len({a['row']for a in terms})==len(terms)>0
        total=[[]for _ in range(4)]
        for a in terms:
            i=a['row'];sign=a['sign'];weight=a['weight']
            assert isinstance(i,int)and 0<=i<12 and sign in (-1,1)
            assert weight and all(isinstance(c,int)for c in weight)
            for j in range(4):total[j]=add(total[j],scale(mul(weight,FORMS[i][j]),sign))
        assert total[:3]==[[],[],[]]
        assert all(isinstance(c,int)for c in cert['residual'])
        assert not add(total[3],cert['residual'])
    applications=0
    for pattern in DATA['patterns']:
        signs=pattern['signature'];assert len(signs)==12 and set(signs)<=set((-1,0,1))
        end=F(0)
        assert pattern['pieces']
        for piece in pattern['pieces']:
            left,right=(F(*x)for x in piece['interval'])
            assert left==end and 0<=left<right<=1
            cert=identities[piece['identity']]
            total=[]
            for a in cert['terms']:
                assert signs[a['row']]==a['sign']
                assert min(bernstein(a['weight'],left,right))>=0
                total=add(total,a['weight'])
            bs=bernstein(total,left,right)
            assert min(bs)>=0 and max(bs)>0
            if left>0:assert value(total,left)>0
            if right<1:assert value(total,right)>0
            assert min(bernstein(cert['residual'],left,right))>=0
            end=right;applications+=1
        assert end==1
    print('EXACTLY VERIFIED',len(identities),'parametric identities and',applications,
        'interval applications covering',len(DATA['patterns']),'sign patterns on (0,1)',flush=True)


if __name__=='__main__':check_certificates()
